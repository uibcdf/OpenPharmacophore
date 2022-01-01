# Open Pharmacophore
from openpharmacophore.io.mol2 import load_mol2_file
from openpharmacophore.utils.random_string import random_string
from openpharmacophore.algorithms.alignment import apply_radii_to_bounds, transform_embeddings
from openpharmacophore._private_tools.exceptions import NoMatchesError, OpenPharmacophoreIOError
from openpharmacophore._private_tools.screening_arguments import check_virtual_screening_kwargs, is_3d_pharmacophore
# Third party
import pandas as pd
from rdkit import RDConfig, Chem, RDLogger, DataStructs
from rdkit.Chem import ChemicalFeatures, rdDistGeom, Descriptors
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D
from rdkit.Chem.Pharm2D.Generate import Gen2DFingerprint
from rdkit.Chem.Pharm3D import EmbedLib
from tqdm.auto import tqdm
# Standard library
import bisect
from collections import namedtuple
import json
from operator import itemgetter
import os

RDLogger.DisableLog('rdApp.*') # Disable rdkit warnings

class VirtualScreening():
    """ Class for performing virtual screening for a database of molecules. 
    
        The database can be fetched, loaded from files or simply passing a list of rdkit molecules.
        Screening can be based on a 3D pharmacophore model or a 2D pharmacophore finegerprint.   

    Parameters
    ----------
    pharmacophore : openpharmacophore.Pharmacophore
        The pharmacophore that will be used to screen the database. Can be a Pharmacophore, 
        StructuredBasedPharmacophore, LigandBasedPharmacophore or a fingerprint

    Attributes
    ----------
    matches : list of 3-tuples (float, str, rdkit.Chem.mol)
        List of molecules that match the pharmacophore. Each tuple is formed by scoring 
        value, the molecule id, and the molecule object.
    
    n_matches : int
        Number of molecules matched to the pharmacophore.

    n_molecules: int
        Number of molcules screened.
    
    n_fails : int
        Number of molecules that cannot be matched to the pharmacophore.
    
    scoring_metric: str
        Metric used to score the molecules, how well they fit to the pharmacophore.

    """
    def __init__(self, pharmacophore, **kwargs):
       
       if is_3d_pharmacophore(pharmacophore):
            self.scoring_metric = "SSD"
            self._screen_fn = self._align_molecules
       elif isinstance(pharmacophore, DataStructs.SparseBitVect): # For pharmacophore fingerprints
            self.scoring_metric = "Similarity"
            self.similarity_fn, self.similarity_cutoff = check_virtual_screening_kwargs(**kwargs)
            self._factory = Gobbi_Pharm2D.factory
            self._screen_fn = self._fingerprint_similarity
       else:
          raise TypeError("pharmacophore must be of type Pharmacophore, StructuredBasedPharmacophore, "
                "LigandBasedPharmacophore, or rdkit.DataStructs.SparseBitVect")

       self.db = ""
       self.matches = []
       self.n_fails = 0
       self.n_matches = 0
       self.n_molecules = 0
       self.pharmacophore = pharmacophore
    
    def get_screening_results(self, form="dataframe"):
        """ Get the results of the screen on a dataframe or a dictionaty

            Parameters
            ---------
            form : {"dataframe", "dict"}
                Whether to return a dataframe or a dictionary.

            Returns
            -------
            pandas.DataFrame or dict
                Dataframe or dictionary with the results information.
        """
        # Values of the scoring that was used for screening. Examples: SSD, 
        # tanimoto similarity
        if self.n_matches == 0:
            raise NoMatchesError("There were no matches in this screen or no database has been screened." 
                                 "Cannot get results.")
    
        score_vals = [i[0] for i in self.matches]
        smiles = [Chem.MolToSmiles(i[2]) for i in self.matches]
        ids = [i[1] for i in self.matches]
        mw = [Descriptors.MolWt(i[2]) for i in self.matches]
        logp = [Descriptors.MolLogP(i[2]) for i in self.matches]

        results = {
            f"{self.db}_id": ids,
            "Smiles": smiles,
            self.scoring_metric: score_vals,
            "Mol_weight": mw,
            "logP": logp
        }

        if form == "dict":
            return results
        elif form == "dataframe":
            return pd.DataFrame().from_dict(results)
        else:
            raise ValueError("form must be dataframe or dict")

    def print_report(self):
        """ Prints a summary report of the screening.
        """
        report_str = self._get_report()
        print(report_str)

    def save_results_to_file(self, file_name):
        """Save the results of the screening to a file. 
        
           The file contains the matched molecules ids, smiles and SSD value. 
           File can be saved as csv or json.

           Parameters
           ----------
           file_name : str
                Name of the file

           Notes
           -----
           Does not return anything. A new file is written.
        """
        file_format = file_name.split(".")[-1]
        if file_format == "csv":
            df = self.get_screening_results(form="dataframe")
            df.to_csv(file_name, index=False)
        elif file_format == "json":
            results = self.get_screening_results(form="dict")
            json_str = json.dumps(results)
            with open(file_name, "w") as f:
                f.write(json_str)
        else:
            raise NotImplementedError
        
    def screen_db_from_dir(self, path, file_extensions=None, **kwargs):
        """ Screen a database of molecules contained in one or more files. 

            Supported file formats can be smi, mol2, sdf.

            Parameters
            ----------
            path : str
                The path to the file or directory that contains the molecules.

            file_extensions : list of str, optional
                A list of file extensions that will be searched for if a directory is passed.
                The default behavior is to load all valid file extensions.

        """
        if not file_extensions:
            file_extensions = ["smi", "mol2", "sdf"]

        if os.path.isdir(path):
            files_list = []
            for root, dirs, files in os.walk(path):
                if '.ipynb_checkpoints' in root:
                    continue
                for file in files:
                    f_extension = file.split(".")[-1]
                    if f_extension not in file_extensions:
                        continue
                    files_list.append(os.path.join(root, file))
            for f in tqdm(files_list):
                molecules = self._load_molecules_file(f, **kwargs)
                self._screen_fn(molecules)

        elif os.path.isfile(path):
            file = path
            molecules = self._load_molecules_file(file, **kwargs)
            self._screen_fn(molecules, pbar=True)

        else:
            raise OpenPharmacophoreIOError("{} is not a valid file/directory".format(path))
    
    def screen_mol_list(self, molecules):
        """Screen a list of molecules

           Parameters
           ----------
           molecules: list of rdkit.Chem.mol
        """
        self._screen_fn(molecules, pbar=True)
    
    def _align_molecules(self, molecules, pbar=False):
        """ Align a list of molecules to a given pharmacophore.

        Parameters
        ----------
        molecules : list of rdkit.Chem.mol
            List of molecules to align.

        verbose : int
            Level of verbosity.
        
        pbar : bool
            If true, a progress bar will be displayed. 
        
        Notes
        -------
        Does not return anything. The attributes matches, n_matches, and n_fails are updated.
        """
        self.n_molecules += len(molecules)

        rdkit_pharmacophore, radii = self.pharmacophore.to_rdkit()
        apply_radii_to_bounds(radii, rdkit_pharmacophore)

        fdef = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
        featFactory = ChemicalFeatures.BuildFeatureFactory(fdef)
        
        Match = namedtuple("Match", ["score", "id", "mol"])
        
        if pbar:
            progress_bar = tqdm(total=self.n_molecules)
        for mol in molecules:

            bounds_matrix = rdDistGeom.GetMoleculeBoundsMatrix(mol)
            # Check if the molecule features can match with the pharmacophore.
            # TODO: replace this function with a custom one that can take other feature definitions
            can_match, all_matches = EmbedLib.MatchPharmacophoreToMol(mol, featFactory, rdkit_pharmacophore)
            # all_matches is a list of tuples where each tuple contains the chemical features
            if can_match:
                # Match the molecule to the pharmacophore without aligning it
                failed, bounds_matrix_matched, matched_mols, match_details = EmbedLib.MatchPharmacophore(all_matches, 
                                                                                                bounds_matrix,
                                                                                                rdkit_pharmacophore, 
                                                                                                useDownsampling=True)
                if failed:
                    self.n_fails += 1
                    if pbar:
                        progress_bar.update()
                    continue
            else:
                self.n_fails += 1
                if pbar:
                    progress_bar.update()
                continue
            atom_match = [list(x.GetAtomIds()) for x in matched_mols]
            try:
                mol_H = Chem.AddHs(mol)
                # Embed molecule onto the pharmacophore
                # embeddings is a list of molecules with a single conformer
                b_matrix, embeddings, num_fail = EmbedLib.EmbedPharmacophore(mol_H, atom_match, rdkit_pharmacophore, count=10)
            except Exception as e:
                self.n_fails += 1
                if pbar:
                    progress_bar.update()
                continue
            # Align embeddings to the pharmacophore 
            SSDs = transform_embeddings(rdkit_pharmacophore, embeddings, atom_match) 
            if len(SSDs) == 0:
                if pbar:
                    progress_bar.update()
                continue
            best_fit_index = min(enumerate(SSDs), key=itemgetter(1))[0]
            
            try:
                mol_id = mol.GetProp("_Name")
            except:
                mol_id = None
            matched_mol = Match(SSDs[best_fit_index], mol_id, embeddings[best_fit_index]) 
            # Append to list in ordered manner
            try:
                bisect.insort(self.matches, matched_mol) 
                self.n_matches += 1
            except:
                # Case when a molecule is repeated. It will throw an error since bisect
                # cannot compare molecules.
                self.n_molecules -= 1
                if pbar:
                    progress_bar.update()
                continue
            if pbar:
                progress_bar.update()

    def _fingerprint_similarity(self, molecules, pbar=False):
        """ Compute fingerprints and similarity values for a list of molecules. 

        Parameters
        ----------
        molecules : list of rdkit.Chem.mol
            List of molecules whose similarity to the pharmacophoric fingerprint will be calculated.
        
        pbar : bool
            If true, a progress bar will be displayed. 
        
        Notes
        -----
        Does not return anything. The attributes matches, n_matches, and n_fails are updated.

        """
        if self.similarity_fn == "tanimoto":
            similarity_fn = DataStructs.TanimotoSimilarity
        elif self.similarity_fn == "dice":
            similarity_fn = DataStructs.DiceSimilarity
       
        Match = namedtuple("Match", ["score", "id", "mol"])
       
        if pbar:
            progress_bar = tqdm(total=len(molecules))
        for mol in molecules:
            self.n_molecules += 1
            fingerprint = Gen2DFingerprint(mol, self._factory)
            similarity = similarity_fn(self.pharmacophore, fingerprint)
            
            if similarity >= self.similarity_cutoff:
                try:
                    mol_id = mol.GetProp("_Name")
                except:
                    mol_id = None
                matched_mol = Match(similarity, mol_id, mol)
                # Append to list in ordered manner
                try:
                    bisect.insort(self.matches, matched_mol)
                    self.n_matches += 1
                except:
                    # Case when a molecule is repeated. It will throw an error since bisect
                    # cannot compare molecules.
                    self.n_molecules -= 1
                    continue
            else:
                self.n_fails += 1 
            
            if pbar:
                progress_bar.update()

    def _get_pharmacophore_fingerprint(self, molecule):
        """ Compute a pharmacophore fingerprint for a single molecule

            Parameters
            ----------
            molecule : rdkit.Chem.mol
                The molecule whose pharmacophore fingerprint will be computed.
            
            Returns
            -------
            fingerprint : rdkit.DataStructs.SparseBitVect
                The pharmacophoric fingerprint.
        """
        factory = Gobbi_Pharm2D.factory
        fingerprint = Gen2DFingerprint(molecule, factory)
        return fingerprint

    def _get_report(self):
        """ Get a report of the screening results.
            
            Returns
            -------
            report_str : str
                The report in string form.
        """
        report_str = "Virtual Screening Results\n"
        report_str += "-------------------------\n"
        report_str += "\nMolecules scanned: " 
        report_str += "{:,}".format(self.n_molecules).rjust(36)
        report_str += "\nMolecules matched to pharmacophore: " 
        report_str += f"{self.n_matches:,}".rjust(19)
        report_str += "\nMolecules that didn't match the pharmacophore: " 
        report_str += "{:,}".format(self.n_fails).rjust(8)
        if self.n_matches > 0:
            report_str += f"\nLowest  {self.scoring_metric} value: "
            report_str += str(round(self.matches[0][0], 4)).rjust(10)
            report_str += f"\nHighest {self.scoring_metric} value: " 
            report_str += str((round(self.matches[-1][0], 4))).rjust(10)
            # Calculate mean
            mean = 0
            N = self.n_matches
            for i in range(N):
                mean += self.matches[i][0]
            mean /= N
            report_str += f"\nAverage {self.scoring_metric} value: " 
            report_str += str(round(mean, 4)).rjust(10)         
            # Print top 5 molecules or less if there are less than 5
            if self.n_matches < 5:
                n_top_mols = min(self.n_matches, 5)
            else:
                n_top_mols = 5                
            report_str += "\n\nTop {} molecules:\n".format(n_top_mols)
            if self.db:
                report_str += f"\n{self.db}ID " + f"{self.scoring_metric}".rjust(12)
            else:
                report_str += "\n   ID   " + f"{self.scoring_metric}".rjust(12)
            report_str += "\n-------".ljust(12) + "------\n".rjust(10)
            for i in range(n_top_mols):
                if self.scoring_metric == "Similarity":
                    i = -(i + 1)
                id = str(self.matches[i][1])
                score = str(round(self.matches[i][0], 4))
                report_str += id.ljust(12)
                report_str += score.rjust(8) + "\n"
        
        return report_str

    def _load_molecules_file(self, file_name, fextension="", **kwargs):
        """ Load a file of molecules into a list of rdkit molecules.

            Parameters
            ----------
            file_name : str
                Name of the file that will be loaded.
            
            fextension : str
                The extension of the file. Must be passed if file argument is a temporary file.
            
            Returns
            -------
            list of rdkit.Chem.mol
        """
        if not fextension:
            fextension = file_name.split(".")[-1]
        
        if fextension == "smi":
            if kwargs:
                if "delimiter" in kwargs and "titleLine" in kwargs:
                    ligands = Chem.SmilesMolSupplier(file_name, delimiter=kwargs["delimiter"], titleLine=kwargs["titleLine"])
                elif "delimiter" in kwargs:
                    ligands = Chem.SmilesMolSupplier(file_name, delimiter=kwargs["delimiter"], titleLine=True)
                elif "titleLine" in kwargs:
                    ligands = Chem.SmilesMolSupplier(file_name, delimiter=' ', titleLine=kwargs["titleLine"])
            else:    
                ligands = Chem.SmilesMolSupplier(file_name, delimiter=' ', titleLine=True)
        elif fextension == "mol2":
            ligands = load_mol2_file(file_name)
        elif fextension == "sdf":
            ligands = Chem.SDMolSupplier(file_name)
        else:
            raise NotImplementedError
        
        # For some strange reason list(ligands) doesn't work if len(ligands) is not called before
        len(ligands)
        ligands = list(ligands)
        ligands = [lig for lig in ligands if lig is not None]
        if len(ligands) == 0:
            raise OpenPharmacophoreIOError("Molecules couldn't be loaded")
        
        return ligands
            
    def __repr__(self):
        return (f"{self.__class__.__name__}(n_matches={self.n_matches}; "
               f"n_fails={self.n_fails})")
        
class ZincMultiScreening():
    """ Class to perform screening of ZINC with multiple pharmacophores.
    
        Parameters
        ----------
        screeners: list of openpharmacophore.VirtualScreening
            A list with virtual screening objects.
        
        download_path: str, default=""
            The path where the files will be saved. If nothing is passed files will
            be deleted afterwards. 
    """
    def __init__(self, screeners, download_path=""):
        if not isinstance(screeners, list):
            raise TypeError("screeners must be of type list.")
        self.screeners = screeners
        if not download_path:
            self.download_path = "./tmp" + random_string(10)
            self.remove_files = True
        else:
            self.download_path = download_path
            self.remove_files = False
        self.download = False
        self.files = []

    def __enter__(self):
        self.download = True
        return self

    def screen(self, subset="Lead-Like", mw_range=None, logp_range=None):
        """ Screen the ZINC subset with each pharmacophore
        
            Parameters
            ----------
            subset : {"Drug-Like", "Lead-Like", "Lugs", "Goldilocks", "Fragments", "Flagments", 
                    "Big-n-Greasy", "Shards"}
                Name of the predifined ZINC subset to be downloaded. 
            
            mw_range : 2-tuple of float, optional 
                Range of molecular weight for the downloaded molecules.
        
            logp_range : 2-tuple of float, optional
                Range of logP for the downloaded molecules.

        """
        for ii, screen in enumerate(self.screeners):
            print(f"Screening with pharmacophore {ii + 1}")
            if self.download:
                self.files = screen.screen_ZINC(
                    download_path=self.download_path,
                    subset=subset,
                    mw_range=mw_range,
                    logp_range=logp_range
                )
                self.download = False
                continue
            screen.screen_db_from_dir(self.download_path)

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.remove_files:
            for file in self.files:
                os.remove(file)
            try:
                os.rmdir(self.download_path)
            except:
                pass
