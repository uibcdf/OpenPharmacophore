from openpharmacophore.databases import chembl
from openpharmacophore.databases.zinc import get_zinc_urls
from openpharmacophore.io.mol2 import load_mol2_file
from openpharmacophore._private_tools.exceptions import MissingParameters
from rdkit import Chem
import pandas as pd
from tqdm.auto import tqdm
import json
import os
from queue import Queue
import requests
import threading
import time

class VirtualScreening():
    """ Base class for performing virtual screening for a database of 
        molecules. The database can be fetched, loaded from a file or
        simply passing a list of rdkit molecules.  

    Parameters
    ----------
    pharmacophore: openpharmacophore.Pharmacophore
        The pharmacophore that will be used to screen the database.

    Attributes
    ----------
    matched_mols: list of rdkit.Chem.mol
        List of molecules that match the pharmacophore.

    n_molecules: int
        Number of molcules screened.

    """

    def __init__(self, pharmacophore):
       self.db = ""
       self.matches = []
       self.scoring_metric = ""
       self.n_fails = 0
       self.n_matches = 0
       self.n_molecules = 0
       self.pharmacophore = pharmacophore
       self._screen_fn = None
       self._file_queue = Queue()
    
    def print_report(self):
        """ Prints a summary report of the screening.
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
            # Calculate mean SSD
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
            report_str += f"\n{self.db}ID " + f"{self.scoring_metric}".rjust(7)
            report_str += "\n-------  " + " ------\n"
            for i in range(n_top_mols):
                report_str += str(self.matches[i][-1]) + "   "
                report_str += str(round(self.matches[i][0], 4)) + "\n"
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
        #TODO: add molecular properties. MW, logP
        file_format = file_name.split(".")[-1]

        # Values of the scoring that was used for screening. Examples: SSD, 
        # tanimoto similarity
        score_vals = [i[0] for i in self.matches]
        smiles = [Chem.MolToSmiles(i[1]) for i in self.matches]
        ids = [i[2] for i in self.matches]

        results = {
            f"{self.db}_id": ids,
            "Smiles": smiles,
            self.scoring_metric: score_vals,
        }

        if file_format == "csv":
            df = pd.DataFrame().from_dict(results)
            df.to_csv(file_name, index=False)
        elif file_format == "json":
            json_str = json.dumps(results)
            with open(file_name, "w") as f:
                f.write(json_str)
        else:
            raise NotImplementedError

    def screen_db(self, db="zinc", download_path=None, **kwargs):
        """ Screen ZINC or ChemBl databases.
            
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
            download_path = "./tmp"
            if not os.path.isdir(download_path):    
                os.mkdir(download_path)
        else:
            delete_files = False

        if db == "zinc":
            self.db = "ZINC"
            if kwargs:
                subset = kwargs["subset"]
                if subset is None:
                    try:
                        mw_range = kwargs["mw_range"]
                        logp_range = kwargs["logp_range"]
                    except:
                        raise MissingParameters("Must pass a molecular weight range and a logP range if no subset is selected") 
                else:
                    mw_range = None
                    logp_range = None
                urls = get_zinc_urls(subset=None, mw_range=mw_range, logp_range=logp_range)
            else:
                urls = get_zinc_urls(subset="Lead-Like")
        
            print("Downloading from ZINC...")
            file_path = self._download_zinc_file(urls[0], download_path)
            self._file_queue.put(file_path)
           
            n_files = len(urls)
            # Start the thread to process files
            threading.Thread(target=self._process_files, daemon=True, args=[delete_files, n_files]).start()

            for url in tqdm(urls[1:]):
                file_path = self._download_zinc_file(url, download_path)
                self._file_queue.put(file_path)
            
            self._file_queue.join()
            print("Finished screening ZINC database")

            try:
                os.rmdir(download_path)
            except:
                pass

        # elif db == "chembl":
        #     self.db = "ChemBL"
        #     # If a subset kwarg is not passed. Use RO5 subset
        #     if kwargs:
        #         subset = kwargs["subset"]
        #     else:
        #         subset = "Ro5"

        #     if subset == "Ro5":
        #         chembl.get_ro5_dataset(download_path)
        #     else:
        #         raise NotImplementedError
        else:
            raise NotImplementedError
        
           
    def screen_db_from_dir(self, path, file_extensions=None):
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
                molecules = self._load_molecules_file(f)
                self._screen_fn(molecules)

        elif os.path.isfile(path):
            file = path
            molecules = self._load_molecules_file(file)
            self._screen_fn(molecules)
            print("File scanned!")

        else:
            raise Exception("{} is not a valid file/directory".format(path))
            
    def _download_zinc_file(self, url, download_path):
        """Download a single file form ZINC.
           
           Parameters
           ----------
           url: str
                URL of the file that will be downloaded.

           download_path: str
                Path where the file will be stored
        """
        try:
            r = requests.get(url, allow_redirects=True)
        except:
            print("Could not download file from {}".format(url))

        file_name =  url[-8:]
        file_path = os.path.join(download_path, file_name)
        with open(file_path, "wb") as file:
            file.write(r.content)
        
        return file_path

    def _download_chembl_file(self, download_path):
        """Download a single file form ChemBl.
           
           Parameters
           ----------

           download_path: str
                Path where the file will be stored
        """
        pass

    def _load_molecules_file(self, file_name):
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
                NUmber of files that will be downloaded

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


class RetrospectiveScreening(VirtualScreening):
    """ Base class for performing retrospective virtual screening. This
        class expects molecules classified as actives and inactives. 

        With this class pharmacophore models can be validated.

    Parameters
    ----------

    Attributes
    ----------

    """

    def __init__(self, pharmacophore):
        self.n_actives = 0
        self.n_inactives = 0
        self.enrichment_factor = 0
        self.auc_score = 0

    def enrichment_plot():
        pass

    def enrichment_factor():
        pass

    def ROC_plot():
        pass

    def AUC():
        pass
