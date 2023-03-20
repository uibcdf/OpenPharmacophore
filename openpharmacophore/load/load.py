from rdkit import Chem

from openpharmacophore import Ligand, Protein
from openpharmacophore.io import load_mol2_ligands, load_sdf
from openpharmacophore.molecular_systems import create_topology
import openpharmacophore.constants as config


class InvalidFileFormatError(ValueError):
    pass


def protein_from_file(traj_file, topology_file):
    """ Create a protein object from a file
    """
    topology, coords = create_topology(traj_file, topology_file)
    return Protein(topology, coords)


def load_ligands_from_file(file_path, file_format):
    """ Load ligands from a file.

        Returns
        -------
        list[Ligand]
    """
    if file_format == "smi":
        return [Ligand(mol) for mol in Chem.SmilesMolSupplier(file_path)]
    if file_format == "sdf":
        return [Ligand(mol) for mol in load_sdf(file_path)]
    if file_format == "mol2":
        return [Ligand(mol) for mol in load_mol2_ligands(file_path)]
    if file_format == "xyz":
        return [Ligand(Chem.MolFromXYZFile(file_path))]
    if file_format == "mol":
        return [Ligand(Chem.MolFromMolFile(file_path))]


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
        return load_ligands_from_file(file_name, file_format)
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
        list[Ligand]
            The set with all the ligands.
    """
    if form not in config.MOL_STR_FORMATS:
        raise InvalidFormError(f"Form {form} is not a supported form")

    return [
        Ligand.from_string(lig, form=form) for lig in ligands
    ]
