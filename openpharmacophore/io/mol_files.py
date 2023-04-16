from rdkit import Chem


# MOL2 files

def mol2(file_name):
    """ Load molecules from a mol2 file.

        Parameters
        ----------
        file_name : str
            Name of the file containing the ligands

        Returns
        ---------
        molecules : list[rdkit.Mol]

    """
    molecules = []
    with open(file_name, 'r') as f:
        doc = [line for line in f.readlines()]

    start = [index for (index, p) in enumerate(doc) if '@<TRIPOS>MOLECULE' in p]
    finish = [index - 1 for (index, p) in enumerate(doc) if '@<TRIPOS>MOLECULE' in p]
    finish.append(len(doc))

    interval = list(zip(start, finish[1:]))

    for ii in interval:
        block = ",".join(doc[ii[0]: ii[1]]).replace(',', '')
        mol = Chem.MolFromMol2Block(block)
        molecules.append(mol)

    return molecules


def iter_mol2(file_name):
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
    with open(file_name) as file_object:
        while True:

            line = file_object.readline()
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
                line = file_object.readline()
                if '@<TRIPOS>MOLECULE' in line or len(line) == 0:
                    break
                mol2_block += line

            yield Chem.MolFromMol2Block(mol2_block)


# SDF files

def sdf(file_path):
    """ Load an sdf file with molecules that may contain multiple conformers.

        Parameters
        ----------
        file_path : str

        Returns
        -------
        list : [rdkit.Chem.Mol]
    """
    supp = Chem.SDMolSupplier(file_path, removeHs=False)
    molecules = {}

    for mol in supp:
        name = mol.GetProp("_Name")
        try:
            molecules[name].AddConformer(mol.GetConformer(), assignId=True)
        except KeyError:
            molecules[name] = mol

    return list(molecules.values())


def iter_sdf(file_):
    raise NotImplementedError


# SMI files

def iter_smi(file_):
    raise NotImplementedError
