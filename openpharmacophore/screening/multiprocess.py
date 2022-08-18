# OpenPharmacophore
from openpharmacophore.screening.screening import VirtualScreening, Match
from openpharmacophore.io.mol_suppliers import smiles_mol_generator, smi_has_header_and_id, mol2_mol_generator
from openpharmacophore._private_tools.exceptions import OpenPharmacophoreIOError, OpenPharmacophoreNotImplementedError
# Third party
from rdkit import Chem
# Standard Library
from multiprocessing import Pool, Queue, Manager, Value
import os
from typing import Callable, List, Tuple


class MultiProcessVirtualScreening(VirtualScreening):
    """ Class to perform virtual screening using multiple cores.
        It generally performs faster than the other
    """

    def __init__(self, pharmacophore, **kwargs):
        super().__init__(pharmacophore, **kwargs)

    @staticmethod
    def _get_files(path: str) -> Queue:
        """ List all files from a directory and put them in a queue.

            Parameters
            ----------
            path : str
                The path of the files.

            Returns
            -------
            file_queue : multiprocessing.Queue
                The file queue.
        """
        valid_formats = ("smi", "txt", "sdf", "mol2", "db2")

        if not os.path.isdir(path):
            raise OpenPharmacophoreIOError("{} is not a valid directory".format(path))

        exclude_prefixes = ('__', '.')

        file_queue = Queue()
        for root, dirs, files in os.walk(path):
            # Ignore hidden folders and files
            files = [f for f in files if not f.startswith(exclude_prefixes)]
            dirs[:] = [d for d in dirs if not d.startswith(exclude_prefixes)]
            for file in files:
                if not file.endswith(valid_formats):
                    continue
                file_queue.put(os.path.join(root, file))

        return file_queue

    def screen_db_from_dir(self, path: str, sort: bool = False) -> None:
        """ Screen a database of molecules contained in one or more files.

            Supported file formats are smi, txt, mol2, sdf.

            Parameters
            ----------
            path : str
                The path to the directory that contains the molecules.

            sort : bool, default=False
                Whether to sort the molecules matched to the pharmacophore by score.

        """
        self._reset()
        self.matches, self.n_molecules = self._multiprocess_screening(path, sort)
        self.n_matches = len(self.matches)
        self.n_fails = self.n_molecules - self.n_matches

    def _multiprocess_screening(self, path: str, sort: bool = False) -> Tuple[List[Match], int]:
        """ Screen a database of molecules using multiple processes.

            Parameters
            ----------
            path : str
                The path to the directory that contains the molecules.

            sort : bool, default=False
                Whether to sort the molecules matched to the pharmacophore by score.

        """
        file_queue = self._get_files(path)
        print("Started Virtual Screening")

        with Manager() as manager:

            matches = manager.list()
            n_molecules = Value("i", 0)

            if self.scoring_metric == "SSD":
                args = (self.rdkit_pharmacophore, matches, self._factory, sort)
                screen_fn = self._align_molecule
            else:
                args = (self.pharmacophore, matches, self._factory, self.similarity_fn,
                        self.similarity_cutoff, sort)
                screen_fn = self._fingerprint_similarity

            pool = Pool(None, self._screen_files, (file_queue, n_molecules, screen_fn, args))
            pool.close()
            pool.join()

            if not sort:
                return list(matches), n_molecules.value

            matches = list(matches)
            matches.sort(key=lambda x: x.score)
            return matches, n_molecules.value

    def _screen_files(self, file_queue: Queue, n_molecules: Value,
                      screen_fn: Callable, args: tuple):
        """ Perform virtual screening to a set of molecules contained in a queue of files.

            Parameters
            ----------
            file_queue : multiprocessing.Queue
                The file queue.

            n_molecules : multiprocessing.Value
                A counter for the total number of molecules screened.

            screen_fn : function
                The function used to screen the molecules. Can be alignment or fingerprint similarity.

            args : tuple
                A tuple with the arguments needed for the screening function.

        """
        while not file_queue.empty():
            file_path = file_queue.get()

            if file_path.endswith(("smi", "txt")):
                has_header, has_id = smi_has_header_and_id(file_path)
                with open(file_path, "r") as fp:
                    for mol in smiles_mol_generator(fp, header=has_header, mol_id=has_id):
                        screen_fn(mol, *args)
                        with n_molecules.get_lock():
                            n_molecules.value += 1

            elif file_path.endswith("mol2"):
                with open(file_path, "r") as fp:
                    for mol in mol2_mol_generator(fp):
                        screen_fn(mol, *args)
                        with n_molecules.get_lock():
                            n_molecules.value += 1

            elif file_path.endswith("sdf"):
                for mol in Chem.SDMolSupplier(file_path):
                    screen_fn(mol, *args)
                    with n_molecules.get_lock():
                        n_molecules.value += 1
            else:
                file_extension = file_path.split(".")[-1]
                raise OpenPharmacophoreNotImplementedError(f"{file_extension} format is currently unsupported.")
