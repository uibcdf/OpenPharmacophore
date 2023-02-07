from openpharmacophore import Ligand, LigandSet, Protein
from openpharmacophore.molecular_systems import create_topology, create_ligand_set
import openpharmacophore.constants as config


class InvalidFileFormatError(ValueError):
    pass


def protein_from_file(traj_file, topology_file):
    """ Create a protein object from a file
    """
    topology, coords = create_topology(traj_file, topology_file)
    return Protein(topology, coords)


def load(file_name, topology_file=None):
    """ Load ligands, protein, protein-ligand complexes, or MD trajectories
        from a file.

        Parameters
        ----------
        file_name : str
            Name of the file containing the data.

        topology_file : str, optional
            File with the topology of the system for a MD trajectory.

        Returns
        -------
        LigandSet or Protein
    """
    file_format = file_name.split(".")[-1]

    if file_format in config.TRAJ_FORMATS:

        if topology_file is not None and topology_file not in config.TOP_FORMATS:
            file_format = topology_file.split(".")[-1]
            raise InvalidFileFormatError(f"File format {file_format} is not supported")

        return protein_from_file(file_name, topology_file)

    elif file_format in config.MOL_FORMATS:
        return create_ligand_set(file_name)
    else:
        raise InvalidFileFormatError(f"File format {file_format} is not supported")


class InvalidFormError(ValueError):
    pass


def load_ligands(ligands, form):
    """ Load ligands from a list of SMILES, SMARTS, Inchi, Mol2 block
        or PDB block.

        Parameters
        ----------
        ligands : list[str]
            List with the ligands

        form : str
            The form of the ligands. Can be "smi", "smarts", "inchi",
            "mol2", "pdb".

        Returns
        -------
        LigandSet
            The set with all the ligands.
    """
    if form not in config.MOL_STR_FORMATS:
        raise InvalidFormError(f"Form {form} is not a supported form")

    lig_set = LigandSet()
    for lig in ligands:
        lig_set.add(Ligand.from_string(lig, form))
    return lig_set
