from openpharmacophore import LigandBasedPharmacophore, LigandReceptorPharmacophore, Pharmacophore
from openpharmacophore.io import (load_json_pharmacophore, load_mol2_pharmacophoric_points,
                                  pharmacophoric_points_from_ph4_file, read_ligandscout)


class InvalidFileFormat(ValueError):
    """ Exception raised when a file format is not supported or
        is incorrect.
    """

    def __init__(self, file_format):

        self.message = "Invalid file format"
        if file_format:
            self.message += f" {file_format}.\n"
        else:
            self.message += ".\n"

        super().__init__(self.message)


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
    pharma = LigandReceptorPharmacophore()
    if file_name.endswith(".json"):
        pharma.add_pharmacophore(
           Pharmacophore(points=load_json_pharmacophore(file_name)[0], ref_struct=0)
        )
    elif file_name.endswith(".mol2"):
        # Loads all pharmacophores contained in the mol2 file. It assumes that each pharmacophore
        # corresponds to a different timestep
        point_lists = load_mol2_pharmacophoric_points(file_name)
        for ii, p_list in enumerate(point_lists):
            pharma.add_pharmacophore(
                Pharmacophore(points=p_list, ref_struct=ii)
            )
    elif file_name.endswith(".pml"):
        pharma.add_pharmacophore(
            Pharmacophore(points=read_ligandscout(file_name), ref_struct=0),
        )
    elif file_name.endswith(".ph4"):
        pharma.add_pharmacophore(
            Pharmacophore(points=pharmacophoric_points_from_ph4_file(file_name),
                          ref_struct=0)
        )
    else:
        raise InvalidFileFormat(file_name.split(".")[-1])
    return pharma


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
        pharmacophore.add_pharmacophore(
            Pharmacophore(load_json_pharmacophore(file_name)[0])
        )
    elif file_name.endswith(".mol2"):
        pharmacophore.add_pharmacophore(
            Pharmacophore(load_mol2_pharmacophoric_points(file_name)[0])
        )
    elif file_name.endswith(".pml"):
        pharmacophore.add_pharmacophore(
            Pharmacophore(read_ligandscout(file_name))
        )
    elif file_name.endswith(".ph4"):
        pharmacophore.add_pharmacophore(
            Pharmacophore(pharmacophoric_points_from_ph4_file(file_name))
        )
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
