from openpharmacophore.databases import chembl, pubchem
from openpharmacophore.databases.zinc import get_zinc_urls
from openpharmacophore.io.mol2 import load_mol2_file
from openpharmacophore.utils.random_string import random_string
from openpharmacophore._private_tools.exceptions import OpenPharmacophoreException
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from tqdm.auto import tqdm
import json
import os
from queue import Queue
import requests
import threading
import time

class VirtualScreening():
    """ Base class for performing virtual screening for a database of 
        molecules. The database can be fetched, loaded from files or
        simply passing a list of rdkit molecules.  

    Parameters
    ----------
    pharmacophore: openpharmacophore.Pharmacophore or list of openpharmacophore.Pharmacophore
        The pharmacophore  or pharmacophores that will be used to screen the database.

    Attributes
    ----------
    matches: list of 3-tuples (float, str, rdkit.Chem.mol)
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

    _screen_fn: function
        The function used for screening.
    
    _file_queue: queue.Queue
        A queue of files that will be screened. Used when downloading from 
        ZINC or ChemBl.

    """

    def __init__(self, pharmacophore):
       self.db = ""
       self.matches = []
       self.scoring_metric = ""
       self.n_fails = 0
       self.n_matches = 0
       self.n_molecules = 0
       self.pharmacophore = pharmacophore
       self.n_pharmacophores = 0
       self._screen_fn = None
       self._file_queue = Queue()
    
    def get_screening_results(self, form="dataframe"):
        """ Get the results of the screen on a dataframe or a
            dictionaty

            Returns
            -------
            pandas.DataFrame
                dataframe with the results information.
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
        else:
            return pd.DataFrame().from_dict(results) 

    def print_report(self):
        """ Prints a summary report of the screening.
        """
        report_str = self._get_report()
        print(report_str)

    def save_results_to_file(self, file_name):
        """Save the results of the screening to a file. The file contains the 
           matched molecules ids, smiles and SSD value. File can be saved as 
           csv or json.

           Parameters
           ----------
           file_name: str
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

    def screen_ZINC(self, db="zinc", download_path=None, subset="Lead-Like", mw_range=None, logp_range=None, **kwargs):
        """ Screen ZINC database.
            
            Parameters
            ---------
            db: str
                Name of the database that will be fetched. Can be "zinc" or "chembl".

            download_path: bool (optional)
                Directory where files will be saved. If None, files will be deleted 
                after processing. Defaults to None

        """
        if not download_path:
            delete_files = True
            download_path = "./tmp" + random_string(10)
            if not os.path.isdir(download_path):    
                os.mkdir(download_path)
        else:
            delete_files = False

        self.db = "ZINC"
        urls = get_zinc_urls(subset=subset, mw_range=mw_range, logp_range=logp_range)
    
        print("Downloading from ZINC...")
        file_path = self._download_zinc_file(urls[0], download_path)
        if file_path is not None:
            self._file_queue.put(file_path)
        
        n_files = len(urls)
        # Start the thread to process files
        threading.Thread(target=self._process_files, daemon=True, args=[delete_files, n_files]).start()

        for url in tqdm(urls[1:]):
            file_path = self._download_zinc_file(url, download_path)
            if file_path is not None:
                self._file_queue.put(file_path)
        
        self._file_queue.join()
        print("Finished screening ZINC database")

        try:
            os.rmdir(download_path)
        except:
            pass
        
    def screen_chembl(self, download_path):
        """ Screen ChemBl database"""
        self.db = "ChemBL"
        pass
        
    def screen_db_from_dir(self, path, file_extensions=None, **kwargs):
        """ Screen a database of molecules contained in one or more files. 
            Format can be smi, mol2, sdf.

            Parmeters
            ---------
            path: str
                The path to the file or directory that contains the molecules.

            file_extensions: list of str (optional)
                A list of file extensions that will be searched for if a directory is passed.
                The default behavior is to load all valid file extensions.

            Notes
            --------
            It does not retur anything. The parameters of the VirtualScreening object are updated accordingly. 
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

    def _download_zinc_file(self, url, download_path):
        """Download a single file form ZINC.
           
           Parameters
           ----------
           url: str
                URL of the file that will be downloaded.

           download_path: str
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

    def _download_chembl_file(self, download_path):
        """Download a single file form ChemBl.
           
           Parameters
           ----------
           download_path: str
                Path where the file will be stored
        """
        pass
    
    def _get_report(self):
        """ Get a report of the screening results.
            
            Returns
            -------
            report_str: str
                The report in string form
        """
        report_str = "Virtual Screening Results\n"
        report_str += "-------------------------\n"
        report_str += "\nMolecules scanned: " 
        report_str += "{:,}".format(self.n_molecules).rjust(36)
        report_str += "\nMolecules matched to pharmacophore: " 
        report_str += str(self.n_matches).rjust(19)
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
        """
            Load a file of molecules of any format and return a list of 
            rdkit molecules.

            Parameters
            ----------
            file_name: str
                Name of the file that will be loaded
            
            Returns
            -------
            A list of rdkit.Chem.mol
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
        """
            For each downloaded file, converts it init a list of rdkit molecules
            and then applies the screenig function to each. 

            Parameters
            ----------
            delete_files: bool
                If true, files will be deleted after processing.

            n_files: int
                Number of files that will be downloaded

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



class RetrospectiveScreening():
    """ Base class for performing retrospective virtual screening. This
        class expects molecules classified as actives and inactives. 

        With this class pharmacophore models can be validated.

    Parameters
    ----------

    Attributes
    ----------

    """
    def __init__(self, pharmacophore):
        self.pharmacophore = pharmacophore
        self.n_molecules = 0
        self.n_actives = 0
        self.n_inactives = 0
        self.n_true_actives = 0 # True positives
        self.n_true_inactives = 0 # True negatives
        self.n_false_actives = 0 # False positives
        self.n_false_inactives = 0 # False negatives
        self.mismatched_actives = []
        self.mismatched_inactives = []
        self.scoring_metric = ""
        self._screen_fn = None

    def from_chembl_target_id(self, target_id, pIC50_threshold=6.3):
        """Retrospective screening from bioactivity data fetched 
           from chembl for the specific target.
           
           Parameters
           ----------
           target_id: str
                ChemBl target id.
           
           pIC50_threshold: float
                The cuttoff value from which a molecule is considered active.
           
           """
        actives, inactives = chembl.get_training_data(target_id, pIC50_threshold)
        
        self.db = "PubChem"
        self.from_training_data(actives, inactives)

    def from_training_data(self, actives, inactives):
        """Retrospective screening from a list of active and inactive molecules
        
            Parameters
            ----------
            actives: 2-tuple
                    The first element is a list of the active compounds ids, and
                    the second elment is a list of smiles for the active compounds
                
            inactives: 2-tuple
                The first element is a list of the inactive compounds ids, and
                the second elment is a list of smiles for the inactive compounds

            Returns
            -------
        
        """
        actives_ids, actives_smiles = actives
        inactives_ids, inactives_smiles = inactives

        self.n_actives = len(actives_ids)
        self.n_inactives = len(inactives_ids)
        self.n_molecules = self.n_actives + self.n_inactives

        # TODO: Add ids to molecule object
        mols_actives = [Chem.MolFromSmiles(smi) for smi in actives_smiles]
        mols_inactives = [Chem.MolFromSmiles(smi) for smi in inactives_smiles]
        
        # Active molecules that were found to be active by the pharmacophore
        matched_actives, self.mismatched_actives = self._screen_fn(mols_actives)
        # Inactive molecules that were found to be active by the pharmacophore
        matched_inactives, self.mismatched_inactives = self._screen_fn(mols_inactives)

        self.n_true_actives = len(matched_actives)
        self.n_true_inactives = len(self.mismatched_inactives)
        self.n_false_actives = len(self.mismatched_actives)
        self.n_false_inactives = len(matched_inactives)
    
    def from_pubchem_bioassay_id(self, bioassay_id):
        """ Retrospective screening from a pubchem bioassay.

            Parameters
            ----------
            bioassay_id: int
                PubChem bioassay id. 
        """
        pubchem_client = pubchem.PubChem()
        actives, inactives = pubchem_client.get_assay_training_data(bioassay_id)
        self.db = "Pubchem"
        self.from_training_data(actives, inactives)

    def from_file(self, file_name):
        pass

    def enrichment_plot(self):
        # To compute enrichment plot and factor we need a list of sorted molecules by a scoring
        # metric i.e SSD values, tanimoto similiarity, etc. We also need to know which molecules
        # are actives and which ones inactives.

        # All molecules need to be assigned a scoring metric else it will not be possible to sort
        # them, which is essential for enrichment calculations.

        # For the ROC plot we need an array with the labels of the molecules, 0 being an inactive
        # molecule and 1 an active one. Moreover, we need an array with the scores. Thats why we
        # need to assign scores to every molecule.

        pass

    def ROC_plot(self):
        pass

    def AUC(self):
        pass

    def enrichment_factor(self):
        n = self.n_true_positives + self.n_true_inactives
        TP = self.n_true_actives
        A = self.n_actives
        N = self.n_inactives
        return (TP / n) / (A / N)

    def sensitivity(self):
        TP = self.n_true_actives
        A = self.n_actives
        return TP / A

    def specificity(self):
        TN = self.n_true_inactives
        FP = self.n_false_actives
        return  TN / (TN + FP)

    def yield_of_actives(self):
        n = self.n_true_positives + self.n_true_inactives
        TP = self.n_true_actives
        return TP / n

    def accuracy(self):
        TP = self.n_true_actives
        TN = self.n_true_inactives
        N = self.n_inactives
        return (TP + TN) / N