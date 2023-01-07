from openpharmacophore.io import MolIO, TrajIO
from openpharmacophore import Ligand, LigandSet, Protein
import openpharmacophore.config as config


class InvalidFileFormatError(ValueError):
    pass


def _load_trajectory(traj_file, topology_file=None):
    """ Load a trajectory from a file and return a Protein.

        Parameters
        ----------
        traj_file : str
            Name of the file containing the data.

        topology_file : str, optional
            File with the topology of the system for a MD trajectory.

        Returns
        -------
        Protein
    """
    if topology_file is not None and topology_file not in config.TOP_FORMATS:
        file_format = topology_file.split(".")[-1]
        raise InvalidFileFormatError(f"File format {file_format} is not supported")

    traj_io = TrajIO(traj_file, topology_file)
    with traj_io:
        traj_io.load_data()
    return traj_io.protein


def _load_ligand_file(file_name):
    """ Load ligands from a file

        Parameters
        ----------
        file_name : str
            Name of the file

        Returns
        -------
        LigandSet

    """
    mol_io = MolIO(file_name)
    with mol_io:
        mol_io.load_data()
    return mol_io.ligands


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
        return _load_trajectory(file_name, topology_file)
    elif file_format in config.MOL_FORMATS:
        return _load_ligand_file(file_name)
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
