from rdkit import Chem


# MOL2 files
def _mol2(file_):
    molecules = []
    doc = [line for line in file_.readlines()]

    start = [index for (index, p) in enumerate(doc) if '@<TRIPOS>MOLECULE' in p]
    finish = [index - 1 for (index, p) in enumerate(doc) if '@<TRIPOS>MOLECULE' in p]
    finish.append(len(doc))

    interval = list(zip(start, finish[1:]))

    for ii in interval:
        block = ",".join(doc[ii[0]: ii[1]]).replace(',', '')
        mol = Chem.MolFromMol2Block(block)
        molecules.append(mol)

    return molecules


def read_mol2(file_name):
    """ Load molecules from a mol2 file.

        Parameters
        ----------
        file_name : str
            Name of the file containing the ligands

        Returns
        ---------
        molecules : list[rdkit.Mol]

    """
    with open(file_name, 'r') as fp:
        molecules = _mol2(fp)
    return molecules


def _iter_mol2(file_):
    line = file_.readline()
    while line:

        if not line:
            break

        if not line.strip():
            continue

        # Skip comments
        if line.startswith("#"):
            continue

        if '@<TRIPOS>MOLECULE' not in line:
            mol2_block = "@<TRIPOS>MOLECULE\n"
            mol2_block += line
        else:
            mol2_block = line
        while True:
            line = file_.readline()
            if '@<TRIPOS>MOLECULE' in line or len(line) == 0:
                break
            mol2_block += line

        yield Chem.MolFromMol2Block(mol2_block)
        line = file_.readline()


def mol2_supplier(file_name):
    """ A molecule generator for mol2 files.

        Parameters
        ----------
        file_name : str
            A file object

        Yields
        ------
        mol : rdkit.Mol
            A molecule.
    """
    with open(file_name) as fp:
        return _iter_mol2(fp)


# SDF files

def _sdf(file_, supplier, remove_hs=False):
    supp = supplier(file_, removeHs=remove_hs)
    molecules = {}

    for mol in supp:
        name = mol.GetProp("_Name")
        try:
            molecules[name].AddConformer(mol.GetConformer(), assignId=True)
        except KeyError:
            molecules[name] = mol

    return list(molecules.values())


def read_sdf(file_path, remove_hs=False):
    """ Load an sdf file with molecules that may contain multiple conformers.

        Parameters
        ----------
        file_path : str

        remove_hs : bool, default=False
            Whether to remove the hydrogens from the molecules.

        Returns
        -------
        list : [rdkit.Chem.Mol]
    """
    return _sdf(file_path, Chem.SDMolSupplier, remove_hs)


def _iter_sdf(file_, supplier, remove_hs=False):
    """ Iterate an sdf file

        Parameters
        ----------
        file_ : str or FileIO
            A path to a file or an stream io such as a BytesIO

        supplier : Callable

        remove_hs : bool

        Returns
        -------
        Iterable
            An iterable of molecules.
    """
    return supplier(file_, removeHs=remove_hs)


# SMI files
def _parse_smi(line, sep):
    fragments = line.split(sep)
    mol = Chem.MolFromSmiles(fragments[0])
    if len(fragments) > 1:
        mol.SetProp("_Name", fragments[1])
    return mol


def _smi(file_, sep=None, header=False):
    start = 1 if header else 0
    lines = file_.readlines()
    molecules = []
    for ii in range(start, len(lines)):
        molecules.append(_parse_smi(lines[ii], sep))
    return molecules


def read_smi(file_name, sep=None, header=False):
    """ Read molecules from a smi file.

        Parameters
        ----------
        file_name : str
            Name or path to the file

        sep : str
            Separator between text in a line of the file. Default
            behavior is to split between whitespace.

        header : bool, default=False
            Whether the file contains a header
    """
    with open(file_name) as fp:
        molecules = _smi(fp, sep, header)
    return molecules


def _iter_smi(file_, sep=None, header=False):
    """ Iterate a smi file.

        Parameters
        ----------
        file_ : FileIO
            A file like object.

        sep : str

        header : bool

        Yields
        ------
        rdkit.Mol
    """
    if header:
        file_.readline()

    line = file_.readline()
    while line:
        yield _parse_smi(line, sep)
        line = file_.readline()
