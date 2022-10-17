

def load_from_file(file_name, pharma_type="ligand-receptor"):
    """ Loads a pharmacophore from a file.

        Parameters
        ----------
        file_name : str
            Path to the file containing the pharmacophores

        pharma_type : str, optional
            Type of pharmacophore, can be "ligand-receptor", "ligand" or
            "receptor".

        Returns
        -------
        Pharmacophore
    """
    if pharma_type == "ligand-receptor":
        pass
    elif pharma_type == "ligand":
        pass
    elif pharma_type == "receptor":
        raise NotImplemented
    else:
        raise ValueError
