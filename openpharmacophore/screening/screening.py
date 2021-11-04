# Open Pharmacophore
from openpharmacophore import Pharmacophore, StructuredBasedPharmacophore, LigandBasedPharmacophore
from openpharmacophore.databases.zinc import get_zinc_urls
from openpharmacophore.io.mol2 import load_mol2_file
from openpharmacophore.utils.random_string import random_string
from openpharmacophore.screening.alignment import apply_radii_to_bounds, transform_embeddings
from openpharmacophore._private_tools.exceptions import OpenPharmacophoreException
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
import json
from operator import itemgetter
import os
from queue import Queue
import requests
import threading
import time

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
       
       if (isinstance(pharmacophore, Pharmacophore) 
        or isinstance(pharmacophore, StructuredBasedPharmacophore)
        or isinstance(pharmacophore, LigandBasedPharmacophore)):
            self.scoring_metric = "SSD"
            self._screen_fn = self._align_molecules

       elif isinstance(pharmacophore, DataStructs.SparseBitVect): # For pharmacophore fingerprints
            self.scoring_metric = "Similarity"

            if kwargs:
                if "similarity" in kwargs:
                    if kwargs["similarity"] != "tanimoto" and kwargs["similarity"] != "dice":
                        raise NotImplementedError
                    self.similiarity_fn = kwargs["similarity"]
                else:
                    self.similiarity_fn = "tanimoto"

                if "sim_cutoff" in kwargs:
                    if kwargs["sim_cutoff"] < 0 and kwargs["sim_cutoff"] > 1:
                        raise ValueError("Similarity cutoff value must lie between 0 and 1")
                    self.similarity_cutoff = kwargs["sim_cutoff"]
                else:
                    self.similarity_cutoff = 0.2
            else:
                self.similiarity_fn = "tanimoto"
                self.similarity_cutoff = 0.2

            self._factory = Gobbi_Pharm2D.factory
            self._screen_fn = self._fingerprint_similarity
       else:
        raise TypeError("pharmacophore must be of type Pharmacophore or StructuredBasedPharmacophore or LigandBasedPharmacophore")

       self.db = ""
       self.matches = []
       self.n_fails = 0
       self.n_matches = 0
       self.n_molecules = 0
       self.pharmacophore = pharmacophore
       self._file_queue = Queue()
    
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
            raise OpenPharmacophoreException("""There were no matches in this screen or no database has been screened. 
            Cannot get results.""")
    
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

    def screen_ZINC(self, subset, mw_range=None, logp_range=None, download_path=""):
        """ Screen a subset of ZINC database. 
        
            The subset can be one of ZINC predefined subsets. Alternatively a subset with a custom molecular 
            weight and logP range can be screened.
            
            Parameters
            ---------
            subset : {"Drug-Like", "Lead-Like", "Lugs", "Goldilocks", "Fragments", "Flagments", 
                    "Big-n-Greasy", "Shards"}
                Name of the predifined ZINC subset to be downloaded. 
            
            mw_range : 2-tuple of float, optional 
                Range of molecular weight for the downloaded molecules.
        
            logp_range : 2-tuple of float, optional
                Range of logP for the downloaded molecules.

            download_path : str, default="" 
                Directory where files will be saved. If None, files will be deleted 
                after processing.
            
            Note
            -----
            If mw_range argument is passed, logp_range must also be passed.

        """
        self.db = "ZINC"
        if not download_path:
            delete_files = True
            download_path = "./tmp" + random_string(10)
        else:
            delete_files = False

        if not os.path.isdir(download_path):    
                os.mkdir(download_path)

        urls = get_zinc_urls(subset=subset, mw_range=mw_range, logp_range=logp_range) 
        n_files = len(urls)

        # Start the thread to process files
        threading.Thread(target=self._process_files, daemon=True, args=[delete_files, n_files]).start()
        files = self.download_zinc_subset(urls, download_path)
        self._file_queue.join()

        if delete_files:
            try:
                os.rmdir(download_path)
            except:
                pass
        else:
            return files
        
    def screen_chembl(self, download_path):
        """ Screen ChemBl database"""
        self.db = "ChemBL"
        pass
        
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
            self._screen_fn(molecules)
            print("File scanned!")

        else:
            raise IOError("{} is not a valid file/directory".format(path))
    
    def screen_mol_list(self, molecules, **kwargs):
        """Screen a list of molecules

           Parameters
           ----------
           molecules: list of rdkit.Chem.mol
        """
        self._screen_fn(molecules, **kwargs)

    def _download_chembl_file(self, download_path):
        """Download a single file form ChemBl.
           
           Parameters
           ----------
           download_path: str
                Path where the file will be stored
        """
        pass

    def _download_zinc_file(self, url, download_path):
        """Download a single file form ZINC.
           
           Parameters
           ----------
           url : str
                URL of the file that will be downloaded.

           download_path : str
                Path where the file will be stored
        """
        res = requests.get(url, allow_redirects=True)

        if res.status_code != requests.codes.ok:
            print("Could not fetch file from {}".format(url))
            return 

        file_name =  url[-8:]
        file_path = os.path.join(download_path, file_name)
        with open(file_path, "wb") as file:
            file.write(res.content)
        
        return file_path

    def download_zinc_subset(self, urls, download_path):
        """ Download a subset from ZINC database.

           Parameters
           ----------
           urls : list of str
                List with the URLs of the files that will be downloaded.

           download_path : str
                Path where the file will be stored
        """
        print("Downloading from ZINC...")
        files = []
        file_path = self._download_zinc_file(urls[0], download_path)
        if file_path is not None:
            self._file_queue.put(file_path) 
            files.append(file_path)

        for url in tqdm(urls[1:]):
            file_path = self._download_zinc_file(url, download_path)
            if file_path is not None:
                self._file_queue.put(file_path)
                files.append(file_path)

        return files
    
    def _align_molecules(self, molecules, verbose=0, sort=True, pbar=False):
        """ Align a list of molecules to a given pharmacophore.

        Parameters
        ----------
        molecules : list of rdkit.Chem.mol
            List of molecules to align.

        verbose : int
            Level of verbosity.

        sort : bool
            If true, the list of matched molecules will be sorted in ascending order.
        
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
        
        if pbar:
            progress_bar = tqdm(total=self.n_molecules)
        for i, mol in enumerate(molecules):

            if verbose == 1 and i % 100 == 0 and i != 0:
                print(f"Screened {i} molecules. Number of matches: {self.n_matches}; Number of fails: {self.n_fails}")

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
                    if verbose == 2:
                        print(f"Couldn't embed molecule {i}")
                    self.n_fails += 1
                    if pbar:
                        progress_bar.update()
                    continue
            else:
                if verbose == 2:
                    print(f"Couldn't match molecule {i}")
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
                if verbose == 2:
                    print(e)
                    print (f"Bounds smoothing failed for molecule {i}")
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
            matched_mol = (SSDs[best_fit_index], mol_id, embeddings[best_fit_index])
            if sort:
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
            else:
                self.matches.append(matched_mols)
            
            if pbar:
                progress_bar.update()

    def _fingerprint_similarity(self, molecules, sort=True, pbar=False):
        """ Compute fingerprints and similarity values for a list of molecules. 

        Parameters
        ----------
        molecules : list of rdkit.Chem.mol
            List of molecules whose similarity to the pharmacophoric fingerprint will be calculated.
        
        sort : bool
            If true, the list of matcehd molecules will be sorted in ascending order.
        
        pbar : bool
            If true, a progress bar will be displayed. 
        
        Notes
        -----
        Does not return anything. The attributes matches, n_matches, and n_fails are updated.

        """
        if self.similiarity_fn == "tanimoto":
            similarity_fn = DataStructs.TanimotoSimilarity
        elif self.similiarity_fn == "dice":
            similarity_fn = DataStructs.DiceSimilarity
        else:
            raise NotImplementedError

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
                matched_mol = (similarity, mol_id, mol)
                if sort:
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
                    self.matches.append(matched_mol)
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

    def _load_molecules_file(self, file_name, **kwargs):
        """ Load a file of molecules into a list of rdkit molecules.

            Parameters
            ----------
            file_name : str
                Name of the file that will be loaded
            
            Returns
            -------
            list of rdkit.Chem.mol
        """
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
            raise Exception("Molecules couldn´t be loaded")
        
        return ligands
    
    def _process_files(self, delete_files, n_files):
        """ Process a file queue. 

            For each downloaded file, converts it into a list of rdkit molecules
            and then applies the screening function to each. 

            Parameters
            ----------
            delete_files : bool
                If true, files will be deleted after processing.

            n_files : int
                Number of files that will be downloaded.

        """
        # Sleep a little so that the download bar appears first.
        time.sleep(2)
        print("Processing files...")
        pbar = tqdm(total=n_files)
        while True:
            file = self._file_queue.get()
            mols = self._load_molecules_file(file)
            self._screen_fn(mols)
            if delete_files:
                os.remove(file)
            pbar.update()
            self._file_queue.task_done()

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
