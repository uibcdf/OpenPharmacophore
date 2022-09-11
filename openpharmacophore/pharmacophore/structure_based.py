from .pharmacophore import Pharmacophore
from .rdkit_pharmacophore import rdkit_pharmacophore
from ..io import (json_pharmacophoric_elements, ligandscout_xml_tree,
                  mol2_file_info, ph4_string)
from ..io import (load_json_pharmacophore, load_mol2_pharmacophoric_points,
                  pharmacophoric_points_from_ph4_file, read_ligandscout)
import nglview as nv
import json


class StructureBasedPharmacophore(Pharmacophore):
    """ Class to store, and extract pharmacophores from protein-ligand complexes.

        The pharmacophores can be extracted from a pdb file or from a molecular
        dynamics simulation.

    """

    def __init__(self, receptor):
        if receptor is None:
            # Pharmacophores will be stored as a list of pharmacophoric points.
            # A list for each pharmacophore
            self._pharmacophores = []
            self._pharmacophores_frames = []  # Contains the frame to which each pharmacophore belongs
            self._topology = None
            self._coordinates = None
            self._num_frames = 0

    @property
    def num_frames(self):
        return self._num_frames

    @classmethod
    def from_file(cls, file_name):
        pass

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
        """ Add a pharmacophoric point to a pharmacophore in an specific frame."""
        self._pharmacophores[frame].append(point)

    def remove_point(self, index, frame):
        """ Removes a pharmacophoric point from the pharmacophore at the given frame."""
        self._pharmacophores[frame].pop(index)

    def remove_picked_point(self, view):
        pass

    def edit_picked_point(self, view):
        pass

    def add_point_in_picked_location(self, view):
        pass

    def add_to_view(self, view, frame=None):
        """ Add pharmacophore(s) to a ngl view.
        """
        if frame is None or isinstance(frame, list):
            raise NotImplementedError
        else:
            for point in self[frame]:
                point.add_to_ngl_view(view)

    def show(self, frame=None):
        """ Shows a 3D representation of the pharmacophore model. """
        view = nv.NGLWidget()
        self.add_to_view(view, frame)

    def to_json(self, file_name, frame=None):
        """ Save pharmacophore(s) to a json file.
        """
        if frame is None or isinstance(frame, list):
            raise NotImplementedError
        else:
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
        if frame is None or isinstance(frame, list):
            raise NotImplementedError
        pharmacophore_data = mol2_file_info([self[frame]])
        with open(file_name, "w") as fp:
            fp.writelines(pharmacophore_data[0])

    def to_rdkit(self, frame):
        """ Transform a pharmacophore at a given frame to a rdkit pharmacophore.
        """
        return rdkit_pharmacophore(self[frame])

    def __len__(self):
        return len(self._pharmacophores)

    def __getitem__(self, frame):
        return self._pharmacophores[frame]
