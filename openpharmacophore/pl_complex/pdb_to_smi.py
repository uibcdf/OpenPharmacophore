import pickle
from pathlib import Path


def pdb_id_mapper():
    """ Loads the dictionary that maps PDB ligand ids to their
        respective smiles.

        Returns
        -------
        dict[str, str]
    """
    path = Path(__file__).parent / "pdb_to_smi.pickle"
    with open(path, "rb") as fp:
        mapper = pickle.load(fp)
    return mapper
