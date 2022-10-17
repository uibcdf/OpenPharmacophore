from .ligand_based import LigandBasedPharmacophore
from .ligand_receptor import LigandReceptorPharmacophore
from ..io import (load_json_pharmacophore, load_mol2_pharmacophoric_points,
                  pharmacophoric_points_from_ph4_file, read_ligandscout)
from .._private_tools.exceptions import InvalidFileFormat, PDBFetchError


def ligand_receptor_from_pharma_file(file_name):
    """ Load a ligand-receptor pharmacophore from a file.

      Parameters
      ---------
      file_name : str
          Name of the file containing the pharmacophore

      Returns
      -------
      pharmacophore : LigandReceptorPharmacophore

    """
    pharmacophore = LigandReceptorPharmacophore()
    if file_name.endswith(".json"):
        pharmacophore._pharmacophores.append(load_json_pharmacophore(file_name)[0])
        pharmacophore._num_frames = 1
        pharmacophore._pharmacophores_frames.append(0)
    elif file_name.endswith(".mol2"):
        # Loads all pharmacophores contained in the mol2 file. It assumes that each pharmacophore
        # corresponds to a different timestep
        pharmacophore._pharmacophores = load_mol2_pharmacophoric_points(file_name)
        pharmacophore._pharmacophores_frames = list(range(0, len(pharmacophore._pharmacophores)))
        pharmacophore._num_frames = len(pharmacophore._pharmacophores)
    elif file_name.endswith(".pml"):
        pharmacophore._pharmacophores.append(read_ligandscout(file_name))
        pharmacophore._num_frames = 1
        pharmacophore._pharmacophores_frames.append(0)
    elif file_name.endswith(".ph4"):
        pharmacophore._pharmacophores.append(pharmacophoric_points_from_ph4_file(file_name))
        pharmacophore._num_frames = 1
        pharmacophore._pharmacophores_frames.append(0)
    else:
        raise InvalidFileFormat(file_name.split(".")[-1])
    return pharmacophore


def ligand_based_from_file(file_name):
    """ Load a ligand-based pharmacophore from a file.

         Parameters
         ---------
         file_name : str
             Name of the file containing the pharmacophore

         Returns
         -------
         pharmacophore : LigandBasedPharmacophore

    """
    pharmacophore = LigandBasedPharmacophore()
    if file_name.endswith(".json"):
        pharmacophore._points = load_json_pharmacophore(file_name)[0]
    elif file_name.endswith(".mol2"):
        pharmacophore._points = load_mol2_pharmacophoric_points(file_name)[0]
    elif file_name.endswith(".pml"):
        pharmacophore._points = read_ligandscout(file_name)
    elif file_name.endswith(".ph4"):
        pharmacophore._points = pharmacophoric_points_from_ph4_file(file_name)
    else:
        raise InvalidFileFormat(file_name.split(".")[-1])
    return pharmacophore


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
        return ligand_receptor_from_pharma_file(file_name)
    elif pharma_type == "ligand":
        return ligand_based_from_file(file_name)
    elif pharma_type == "receptor":
        raise NotImplemented
    else:
        raise ValueError
