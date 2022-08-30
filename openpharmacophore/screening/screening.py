# OpenPharmacophore
from openpharmacophore import Pharmacophore, LigandBasedPharmacophore, StructuredBasedPharmacophore
from openpharmacophore.io.mol_suppliers import smiles_mol_generator, smi_has_header_and_id, mol2_mol_generator
from openpharmacophore.screening.alignment import apply_radii_to_bounds, transform_embeddings
from openpharmacophore._private_tools.exceptions import NoMatchesError, OpenPharmacophoreIOError, \
    OpenPharmacophoreNotImplementedError, OpenPharmacophoreTypeError
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
from typing import Callable, List, Tuple, TypeVar

RDLogger.DisableLog('rdApp.*')  # Disable rdkit warnings

Match = namedtuple("Match", ["score", "id", "mol"])
PharmacophoreType = TypeVar(
    "PharmacophoreType", LigandBasedPharmacophore,
    StructuredBasedPharmacophore, Pharmacophore, DataStructs.SparseBitVect
)
Results = TypeVar("Results", pd.DataFrame, dict)


class VirtualScreening:
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
    matches : list of 3-tuples (float, str, rdkit.Chem.Mol)
        List of molecules that match the pharmacophore. Each tuple is formed by scoring 
        value, the molecule id, and the molecule object.
    
    n_matches : int
        Number of molecules matched to the pharmacophore.

    n_molecules: int
        Number of molcules screened.
    
    n_fails : int
        Number of molecules that cannot be matched to the pharmacophore.
    
    scoring_metric : str
        Metric used to score the molecules, how well they fit to the pharmacophore.
    
    pharmacophore : openpharmacophore.Pharmacophore
        The pharmacophore that will be used to screen the database. Can be a Pharmacophore, 
        StructuredBasedPharmacophore, LigandBasedPharmacophore or a fingerprint.

    """

    def __init__(self, pharmacophore: PharmacophoreType, **kwargs) -> None:

        if is_3d_pharmacophore(pharmacophore):
            self.scoring_metric = "SSD"
            self._screen_fn = self._align_molecule
            self._factory = ChemicalFeatures.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir,
                                                                              'BaseFeatures.fdef'))
            self.rdkit_pharmacophore = self._get_rdkit_pharmacophore(pharmacophore)
        elif isinstance(pharmacophore, DataStructs.SparseBitVect):  # For pharmacophore fingerprints
            self.scoring_metric = "Similarity"
            self.similarity_fn, self.similarity_cutoff = check_virtual_screening_kwargs(**kwargs)
            self._factory = Gobbi_Pharm2D.factory
            self._screen_fn = self._fingerprint_similarity
        else:
            raise OpenPharmacophoreTypeError(
                "pharmacophore must be of type Pharmacophore, StructuredBasedPharmacophore, "
                "LigandBasedPharmacophore, or rdkit.DataStructs.SparseBitVect")

        self.matches = []
        self.n_fails = 0
        self.n_matches = 0
        self.n_molecules = 0
        self.pharmacophore = pharmacophore

    def get_screening_results(self, form: str = "dataframe") -> Results:
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
            f"Id": ids,
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

    def print_report(self) -> None:
        """ Prints a summary report of the screening.
        """
        report_str = self._get_report()
        print(report_str)

    def save_results_to_file(self, file_name: str) -> None:
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

    def screen_db_from_dir(self, path: str, pbar: bool = True) -> None:
        """ Screen a database of molecules contained in one or more files. 

            Supported file formats are smi, txt, mol2, sdf.

            Parameters
            ----------
            path : str
                The path to the directory that contains the molecules.
                
            pbar : bool
                Whether to display a progress bar

        """
        self._reset()
        valid_formats = ("smi", "txt", "sdf", "mol2", "db2")

        if not os.path.isdir(path):
            raise OpenPharmacophoreIOError("{} is not a valid directory".format(path))

        exclude_prefixes = ('__', '.')
        file_list = []
        for root, dirs, files in os.walk(path):
            # Ignore hidden folders and files
            files = [f for f in files if not f.startswith(exclude_prefixes)]
            dirs[:] = [d for d in dirs if not d.startswith(exclude_prefixes)]
            for file in files:
                if not file.endswith(valid_formats):
                    continue
                file_list.append(os.path.join(root, file))

        for file_path in tqdm(file_list, disable=not pbar):
            self.screen_mol_file(file_path)

    def screen_mol_list(self, molecules: List[Chem.Mol], pbar: bool = True) -> None:
        """ Screen a list of molecules

            Parameters
            ----------
            molecules: list of rdkit.Chem.Mol
            
            pbar : bool
                Whether to display a progress bar
        """
        self._reset()
        self.n_molecules = len(molecules)

        if self.scoring_metric == "SSD":
            for mol in tqdm(molecules, disable=not pbar):
                self._align_molecule(mol, self.rdkit_pharmacophore, self.matches, self._factory, sort=True)
        else:
            for mol in tqdm(molecules, disable=not pbar):
                self._fingerprint_similarity(mol, self.pharmacophore, self.matches, self._factory,
                                             self.similarity_fn, self.similarity_cutoff, sort=True)

        self.n_matches = len(self.matches)
        self.n_fails = self.n_molecules - self.n_matches

    def screen_mol_file(self, file_path: str, sort: bool = True) -> None:
        """ Perform virtual screening to a file of molecules.
        
            Parameters
            ----------
            file_path : str
                The path of the file that contains the molecules.
                
            sort : bool, default=True
                Whether to sort the matched molecules by score.
        """
        if self.scoring_metric == "SSD":
            args = (self.rdkit_pharmacophore, self.matches, self._factory, sort)
            screen_fn = self._align_molecule
        else:
            args = (self.pharmacophore, self.matches, self._factory, self.similarity_fn,
                    self.similarity_cutoff, sort)
            screen_fn = self._fingerprint_similarity

        self._screen_file(file_path, screen_fn, args)
        self.n_matches = len(self.matches)
        self.n_fails = self.n_molecules - self.n_matches

    def _screen_file(self, file_path: str, screen_fn: Callable, args: Tuple) -> None:
        """ Apply a screening function to each molecule in a file.
        
            Parameters
            ----------
            file_path : str
                The path of the file.
                
            screen_fn : function
                The function used to screen the molecules. Can be alignment or fingerprint similarity.
                
            args : tuple
                A tuple with the arguments needed for the screening function.
        """

        if file_path.endswith(("smi", "txt")):
            has_header, has_id = smi_has_header_and_id(file_path)
            with open(file_path, "r") as fp:
                for mol in smiles_mol_generator(fp, header=has_header, mol_id=has_id):
                    screen_fn(mol, *args)
                    self.n_molecules += 1
        elif file_path.endswith("mol2"):
            with open(file_path, "r") as fp:
                for mol in mol2_mol_generator(fp):
                    screen_fn(mol, *args)
                    self.n_molecules += 1
        elif file_path.endswith("sdf"):
            for mol in Chem.SDMolSupplier(file_path):
                screen_fn(mol, *args)
                self.n_molecules += 1
        else:
            file_extension = file_path.split(".")[-1]
            raise OpenPharmacophoreNotImplementedError(f"{file_extension} format is currently unsupported.")

    @staticmethod
    def _align_molecule(mol, pharmacophore, matches, featFactory, sort=False):
        """ Align a molecule to a given pharmacophore.
        
            Uses rdkit alignment algorithm

            Parameters
            ----------
            mol : rdkit.Chem.Mol
                Molecule to align.
                
            matches : list
                If a moleculed is matched to the pharmacophore it will be appended to this list.
                
            pharmacophore : rdkit.Chem.Pharm3D.Pharmacophore
                An rdkit pharmacophore

            featFactory : rdkit.Chem.rdMolChemicalFeatures.MolChemicalFeatureFactory
                The feature factory.
            
            sort : bool, default=False
                Whether to sort the list with the matches

        """
        bounds_matrix = rdDistGeom.GetMoleculeBoundsMatrix(mol)
        # Check if the molecule features can match with the pharmacophore.
        can_match, all_matches = EmbedLib.MatchPharmacophoreToMol(mol, featFactory, pharmacophore)
        # all_matches is a list of tuples where each tuple contains the chemical features
        if can_match:
            # Match the molecule to the pharmacophore without aligning it
            failed, bounds_matrix_matched, matched_mols, match_details = EmbedLib.MatchPharmacophore(
                all_matches, bounds_matrix, pharmacophore, useDownsampling=True)
            if failed:
                return
        else:
            return
        atom_match = [list(x.GetAtomIds()) for x in matched_mols]
        try:
            mol_H = Chem.AddHs(mol)
            # Embed molecule onto the pharmacophore
            # embeddings is a list of molecules with a single conformer
            b_matrix, embeddings, num_fail = EmbedLib.EmbedPharmacophore(mol_H, atom_match, pharmacophore, count=10)
        except:
            return
        # Align embeddings to the pharmacophore 
        SSDs = transform_embeddings(pharmacophore, embeddings, atom_match)
        if len(SSDs) == 0:
            return
        best_fit_index = min(enumerate(SSDs), key=itemgetter(1))[0]
        try:
            mol_id = mol.GetProp("_Name")
        except:
            mol_id = None

        matched_mol = Match(SSDs[best_fit_index], mol_id, embeddings[best_fit_index])
        if sort:
            # Append to list in ordered manner
            try:
                # Case when a molecule is repeated. It will throw an error since bisect
                # cannot compare molecules.
                bisect.insort(matches, matched_mol)
            except:
                return
        else:
            matches.append(matched_mol)

    @staticmethod
    def _fingerprint_similarity(mol, pharmacophore_fp, matches, factory, similarity_fn, sim_cutoff, sort=False):
        """ Compute fingerprints similarity values for a molecule. 

        Parameters
        ----------
        mol : rdkit.Chem.Mol
           Molecule whose similarity to the pharmacophoric fingerprint will be calculated.
           
        matches : list
            If a moleculed is matched to the pharmacophore it will be appended to this list.
                
        pharmacophore_fp : rdkit.DataStructs.SparseBitVect
            An pharmacophore fingerprint.

        factory : 
        
        similarity_fn : function
            The similarity function that will be applied to the molecule. Can be tanimoto or dice.
        
        sim_cutoff : float between 0 and 1
            The cutoff value from which a molecule will be considered a match.
        
        sort : bool, default=False
            Whether to sort the list with the matches.
        
        Notes
        -----
        Does not return anything. The attributes matches, n_matches, and n_fails are updated.

        """
        fingerprint = Gen2DFingerprint(mol, factory)
        similarity = similarity_fn(pharmacophore_fp, fingerprint)

        if similarity >= sim_cutoff:
            try:
                mol_id = mol.GetProp("_Name")
            except:
                mol_id = None
            matched_mol = Match(similarity, mol_id, mol)
            if sort:
                # Append to list in ordered manner
                try:
                    # Case when a molecule is repeated. It will throw an error since bisect
                    # cannot compare molecules.
                    bisect.insort(matches, matched_mol)
                except:
                    return
            else:
                matches.append(matched_mol)

    def _get_pharmacophore_fingerprint(self, molecule):
        """ Compute a pharmacophore fingerprint for a single molecule

            Parameters
            ----------
            molecule : rdkit.Chem.Mol
                The molecule whose pharmacophore fingerprint will be computed.
            
            Returns
            -------
            fingerprint : rdkit.DataStructs.SparseBitVect
                The pharmacophoric fingerprint.
        """
        factory = Gobbi_Pharm2D.factory
        fingerprint = Gen2DFingerprint(molecule, factory)
        return fingerprint

    def _get_rdkit_pharmacophore(self, pharmacophore):
        """ Transfor a pharmacophore to an rdkit pharmacophore
        
            Parameters
            ----------
            pharmacophore : openpharmacophore.Pharmacophore
                The pharmacophore.
                
            Returns
            -------
            rdkit.Chem.Pharm3D.Pharmacophore
                An rdkit pharmacophore with the radii of the points applied.
        """
        rdkit_pharmacophore, radii = pharmacophore.to_rdkit()
        apply_radii_to_bounds(radii, rdkit_pharmacophore)

        return rdkit_pharmacophore

    def _get_report(self) -> None:
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
            report_str += "\n   ID   " + f"{self.scoring_metric}".rjust(12)
            report_str += "\n-------".ljust(12) + "------\n".rjust(10)
            for i in range(n_top_mols):
                if self.scoring_metric == "Similarity":
                    i = -(i + 1)
                id_ = str(self.matches[i][1])
                score = str(round(self.matches[i][0], 4))
                report_str += id_.ljust(12)
                report_str += score.rjust(8) + "\n"

        return report_str

    def _reset(self) -> None:
        """ Resets the attributes of the VirtualScreening instance."""
        self.matches.clear()
        self.n_fails = 0
        self.n_matches = 0
        self.n_molecules = 0

    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}(n_matches={self.n_matches}; "
                f"n_fails={self.n_fails})")
