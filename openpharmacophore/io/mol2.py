from rdkit import Chem


def load_mol2_ligands(file_name):
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
