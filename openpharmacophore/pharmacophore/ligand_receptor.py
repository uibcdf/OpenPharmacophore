from .pharmacophore import Pharmacophore
from .pharmacophoric_point import PharmacophoricPoint
from .rdkit_pharmacophore import rdkit_pharmacophore
from ..io import (json_pharmacophoric_elements, ligandscout_xml_tree,
                  mol2_file_info, ph4_string)
from ..io import (load_json_pharmacophore, load_mol2_pharmacophoric_points,
                  pharmacophoric_points_from_ph4_file, read_ligandscout)
from .._private_tools.exceptions import InvalidFileFormat, PDBFetchError
import numpy as np
import nglview as nv
import pyunitwizard as puw
import json
import re
import requests


class LigandReceptorPharmacophore(Pharmacophore):
    """ Class to store, and extract pharmacophores from protein-ligand complexes.

        The pharmacophores can be extracted from a pdb file or from a molecular
        dynamics simulation.

        Parameters
        ----------
        receptor : Any
            A PDB id, pdb file or a trajectory or topology file containing the protein-ligand
            complex.

    """

    def __init__(self, receptor):
        # Pharmacophores will be stored as a list of pharmacophoric points.
        # A list for each pharmacophore
        self._pharmacophores = []
        self._pharmacophores_frames = []  # Contains the frame to which each pharmacophore belongs
        self._num_frames = 0

        self._pdb = ""  # A file path to a pdb file or a pdb as a string

        if isinstance(receptor, str):
            if self._is_pdb_id(receptor):
                self._pdb = self._fetch_pdb(receptor)
                self._num_frames += 1
            elif receptor.endswith(".pdb"):
                self._pdb_is_str = False
                self._pdb = receptor
                self._num_frames += 1
            else:
                # Load from a trajectory
                raise NotImplementedError

    @property
    def num_frames(self):
        return self._num_frames

    @staticmethod
    def _is_pdb_id(receptor):
        """ Check if the receptor is a PDB id.

            Parameters
            ----------
            receptor: str
                The receptor should be a string

            Returns
            -------
            bool
                Whether the receptor is a pdb id.
        """
        if len(receptor) == 4:
            pattern = re.compile('[0-9][a-zA-Z_0-9]{3}')
            if pattern.match(receptor):
                return True
        return False

    @staticmethod
    def _fetch_pdb(pdb_id):
        """ Fetch a PDB with the given id.
        """
        url = f'http://files.rcsb.org/download/{pdb_id}.pdb'
        res = requests.get(url, allow_redirects=True)

        if res.status_code != 200:
            raise PDBFetchError(pdb_id, url)

        return res.content.decode()

    def from_file(self, file_name):
        """ Load a pharmacophore from a file.

          Parameters
          ---------
          file_name : str
              Name of the file containing the pharmacophore

        """
        if file_name.endswith(".json"):
            self._pharmacophores.append(load_json_pharmacophore(file_name)[0])
            self._num_frames = 1
            self._pharmacophores_frames.append(0)
        elif file_name.endswith(".mol2"):
            # Loads all pharmacophores contained in the mol2 file. It assumes that each pharmacophore
            # corresponds to a different timestep
            self._pharmacophores = load_mol2_pharmacophoric_points(file_name)
            self._pharmacophores_frames = list(range(0, len(self._pharmacophores)))
            self._num_frames = len(self._pharmacophores)
        elif file_name.endswith(".pml"):
            self._pharmacophores.append(read_ligandscout(file_name))
            self._num_frames = 1
            self._pharmacophores_frames.append(0)
        elif file_name.endswith(".ph4"):
            self._pharmacophores.append(pharmacophoric_points_from_ph4_file(file_name))
            self._num_frames = 1
            self._pharmacophores_frames.append(0)
        else:
            raise InvalidFileFormat(file_name.split(".")[-1])

    def add_frame(self):
        """ Add a new frame to the pharmacophore. """
        self._pharmacophores.append([])
        self._pharmacophores_frames.append(self._num_frames)
        self._num_frames += 1

    def add_points_to_frame(self, point_list, frame):
        """ Add pharmacophoric points from a list to a frame. """
        for point in point_list:
            self._pharmacophores[frame].append(point)

    def add_point(self, point, frame):
        """ Add a pharmacophoric point to a pharmacophore in a specific frame."""
        self._pharmacophores[frame].append(point)

    def remove_point(self, index, frame):
        """ Removes a pharmacophoric point from the pharmacophore at the given frame."""
        self._pharmacophores[frame].pop(index)

    def remove_picked_point(self, view):
        raise NotImplementedError

    def edit_picked_point(self, view):
        raise NotImplementedError

    def add_point_in_picked_location(self, view):
        raise NotImplementedError

    def add_to_view(self, view, frame=0):
        """ Add pharmacophore(s) to a ngl view.
        """
        if isinstance(frame, list):
            raise NotImplementedError
        else:
            for point in self[frame]:
                point.add_to_ngl_view(view)

    def show(self, frame=0):
        """ Shows a 3D representation of the pharmacophore model. """
        view = nv.NGLWidget()
        self.add_to_view(view, frame)
        return view

    def to_json(self, file_name, frame):
        """ Save pharmacophore(s) to a json file.
        """
        data = json_pharmacophoric_elements(self[frame])
        with open(file_name, "w") as fp:
            json.dump(data, fp)

    def to_ligand_scout(self, file_name, frame):
        """ Save a pharmacophore at a given frame to ligand scout format (pml).
        """
        xml_tree = ligandscout_xml_tree(self[frame])
        xml_tree.write(file_name, encoding="UTF-8", xml_declaration=True)

    def to_moe(self, file_name, frame):
        """ Save a pharmacophore at a given frame to moe format (ph4).
        """
        pharmacophore_str = ph4_string(self[frame])
        with open(file_name, "w") as fp:
            fp.write(pharmacophore_str)

    def to_mol2(self, file_name, frame=None):
        """ Save pharmacophore(s) to mol2 file.
        """
        # TODO: save multiple pharmacophores
        if frame is None or isinstance(frame, list):
            raise NotImplementedError
        pharmacophore_data = mol2_file_info([self[frame]])
        with open(file_name, "w") as fp:
            fp.writelines(pharmacophore_data[0])

    def to_rdkit(self, frame):
        """ Transform a pharmacophore at a given frame to a rdkit pharmacophore.
        """
        return rdkit_pharmacophore(self[frame])

    def extract(self, ligand_id, frames=None):
        """ Extract pharmacophore(s) from the receptor. A protein-ligand complex
            can contain multiple ligands or small molecules, pharmacophore(s) is
            extracted only for the selected one.

            Parameters
            ----------
            ligand_id : str
                The id of the ligand whose pharmacophore will be extracted.

            frames : list[int] or 'all', optional
                Extract pharmacophores from the given frames of the trajectory.
                If None is passed only the first frame will be used.
        """
        raise NotImplementedError

    def __len__(self):
        return len(self._pharmacophores)

    def __getitem__(self, frame):
        return self._pharmacophores[frame]
