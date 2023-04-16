from openpharmacophore import Ligand
from openpharmacophore.io import mol_files, exceptions


class MolDB:
    """ A repository of molecules that can be easily iterated.
    """

    _mol_file_iterators = {
        "smi": mol_files._iter_smi,
        "sdf": mol_files._iter_sdf,
        "mol2": mol_files._iter_mol2,
    }

    def __init__(self):
        self._read_fn = None
        self._molecules = None

    def from_smiles(self, smiles):
        """ Adds molecules from a list of SMILES.

            Parameters
            ----------
            smiles : list[str]
                List with the SMILES.
        """
        self._read_fn = lambda m: Ligand.from_string(m, form="smi")
        self._molecules = smiles

    def from_rdkit(self, mols):
        """ Adds molecules from a list of rdkit Mol.

           Parameters
           ----------
           mols : list[rdkit.Mol]
               List with the molecules.
       """
        self._read_fn = lambda m: Ligand(m)
        self._molecules = mols

    def from_file(self, file_name):
        """ Adds molecules from a file.

            Parameters
            ----------
            file_name: str
                Name or path of the file.
        """
        file_extension = file_name.split(".")[-1]
        with open(file_name, "r") as fp:
            self._from_file(fp, file_extension)

    def _from_file(self, file_, extension):
        self._read_fn = lambda m: Ligand(m)
        try:
            self._molecules = MolDB._mol_file_iterators[extension](file_)
        except KeyError:
            raise exceptions.InvalidFileFormatError(f"{extension} is not a supported file format")

    def _iter_list(self):
        """ Iterate a list of molecules.
        """
        if self._ii >= len(self._molecules):
            del self._ii
            raise StopIteration
        mol = self._read_fn(self._molecules[self._ii])
        self._ii += 1
        return mol

    def __iter__(self):
        if isinstance(self._molecules, (list, tuple)):
            self._ii = 0
        return self

    def __next__(self):
        """ Iterate the molecules of the database converting
            them to the Ligand class.
        """
        if isinstance(self._molecules, (list, tuple)):
            return self._iter_list()
        return next(self._molecules)
