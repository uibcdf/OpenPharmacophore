from openpharmacophore import Ligand
from openpharmacophore.io import mol_files, exceptions


class InMemoryMolDB:
    """ A repository of molecules that are stored in memory. Such as
        a list of smiles or rdkit molecules.
    """

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

    def __iter__(self):
        self._ii = 0
        return self

    def __next__(self):
        """ Iterate the molecules of the database converting
            them to the Ligand class.
        """
        if self._ii >= len(self._molecules):
            del self._ii
            raise StopIteration
        mol = self._read_fn(self._molecules[self._ii])
        self._ii += 1
        return mol


class LigandGenerator:

    def __init__(self, generator):
        self._gen = generator

    def __iter__(self):
        return self

    def __next__(self):
        return Ligand(next(self._gen))


class MolDB:
    """ A repository of molecules that are stored in files or file
        like objects.
    """

    _mol_file_iterators = {
        "smi": mol_files._iter_smi,
        "sdf": mol_files._iter_sdf,
        "mol2": mol_files._iter_mol2,
    }

    def __init__(self):
        self._files = []
        self._extensions = []

    def from_file(self, file_name):
        """ Adds molecules from a file.

            Parameters
            ----------
            file_name: str
                Name or path of the file.
        """
        file_extension = file_name.split(".")[-1]
        if file_extension not in MolDB._mol_file_iterators:
            raise exceptions.InvalidFileFormatError(
                f"{file_extension} is not a supported file format"
            )
        self._from_file(open(file_name), file_extension)

    def from_file_list(self, files):
        """ Adds molecules from a list of files.

            Parameters
            ----------
            files: list[str]
                Name or path of the file.
        """
        for filename in files:
            self.from_file(filename)

    def _from_file(self, file_, extension):
        self._files.append(file_)
        self._extensions.append(extension)

    def __iter__(self):
        for ii in range(len(self._files)):
            file_ = self._files[ii]
            extension = self._extensions[ii]
            yield from LigandGenerator(MolDB._mol_file_iterators[extension](file_))
            file_.close()
