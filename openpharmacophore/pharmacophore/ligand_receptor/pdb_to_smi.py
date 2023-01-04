import pickle


def pdb_id_mapper():
    """ Loads the dictionary that maps PDB ligand ids to their
        respective smiles.

        Returns
        -------
        dict[str, str]
    """
    with open("pdb_to_smi.pickle") as fp:
        mapper = pickle.load(fp)
    return mapper
