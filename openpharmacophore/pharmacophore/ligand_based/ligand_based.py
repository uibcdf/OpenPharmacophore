from openpharmacophore import distance_between_pharmacophoric_points, PharmacophoricPoint
from openpharmacophore.pharmacophore.ligand_based.rdp import find_common_pharmacophores
from openpharmacophore.pharmacophore.pharmacophore import Pharmacophore
from openpharmacophore.pharmacophore.rdkit_pharmacophore import rdkit_pharmacophore
import openpharmacophore.io as io
from openpharmacophore._private_tools.exceptions import InvalidFileFormat
import numpy as np
import nglview as nv
import pyunitwizard as puw
import rdkit.Chem.AllChem as Chem
import json
import os


class LigandBasedPharmacophore(Pharmacophore):
    """ Class to store and extract pharmacophores from a set of ligands.

    """

    def __init__(self):
        self._points = []
        self._ligands = []

    @property
    def pharmacophoric_points(self):
        return self._points

    @pharmacophoric_points.setter
    def pharmacophoric_points(self, points):
        self._points = points

    @property
    def num_points(self):
        return len(self._points)

    @property
    def ligands(self):
        return self._ligands

    @ligands.setter
    def ligands(self, ligand_list):
        self._ligands = ligand_list

    @ligands.deleter
    def ligands(self):
        self._ligands.clear()

    def load_ligands(self, ligands):
        """ Load ligands from a file and store them as a
            list of rdkit molecules.

            Parameters
            -----------
            ligands : str
                A path to a file.
        """
        if isinstance(ligands, str) and os.path.isfile(ligands):
            self._ligands = io.mol_file_to_list(ligands)

    def load_ligands_from_smi(self, ligands):
        """ Load ligands from a list of smiles and store them as a
            list of rdkit molecules.

            Parameters
            -----------
            ligands : list[str]
                A path to a file.
        """
        self._ligands = [Chem.MolFromSmiles(mol) for mol in ligands]

    def from_file(self, file_name):
        """ Load a pharmacophore from a file.

           Parameters
           ---------
           file_name : str
               Name of the file containing the pharmacophore

       """
        if file_name.endswith(".json"):
            self._points = io.load_json_pharmacophore(file_name)[0]
        elif file_name.endswith(".mol2"):
            self._points = io.load_mol2_pharmacophoric_points(file_name)[0]
        elif file_name.endswith(".pml"):
            self._points = io.read_ligandscout(file_name)
        elif file_name.endswith(".ph4"):
            self._points = io.pharmacophoric_points_from_ph4_file(file_name)
        else:
            raise InvalidFileFormat(file_name.split(".")[-1])

    def add_point(self, point):
        """ Adds a pharmacophoric point.

            Parameters
            ----------
            point : PharmacophoricPoint
                The pharmacophoric point that will be added.
        """
        self._points.append(point)

    def remove_point(self, index):
        """ Removes a pharmacophoric point from the pharmacophore.

            Parameters
            ----------
            index: int
                The index of the pharmacophoric point.
        """
        self._points.pop(index)

    @staticmethod
    def get_picked_point_index(view):
        """ Get the index of a point picked in a view in the
            pharmacophore.

            Parameters
            ----------
            view : nglview.NGLWidget
                View where the pharmacophore will be added.

            Returns
            -------
            index: int or None
                The index or None if something invalid was selected.
        """
        if len(view.picked) != 1:
            # An atom or nothing was selected
            return

        picked_index = view.picked["component"]
        # We'll assume that the view components are the molecule(s) followed by
        # the pharmacophoric points
        points_start = None
        for ii in range(len(view._ngl_component_names)):
            if view._ngl_component_names[ii] == "nglview.shape.Shape":
                points_start = ii
                break

        if points_start is None:
            return

        return picked_index - points_start

    def remove_picked_point(self, view):
        """ Remove a pharmacophoric point selected in a ngl view.

            If nothing is selected or the selected thing is not a
            pharmacophoric point, it doesn't change the view.

           Parameters
           ----------
           view : nglview.NGLWidget
               View where the pharmacophore will be added.
        """
        index = self.get_picked_point_index(view)
        if index is not None and index < len(self):
            self.remove_point(index)

    def edit_picked_point(self, view, center, radius):
        """ Remove a pharmacophoric point selected in a ngl view.

           If nothing is selected or the selected thing is not a
           pharmacophoric point, it doesn't change the view.

          Parameters
          ----------
          view : nglview.NGLWidget
              View where the pharmacophore will be added.

          center: Quantity
              The new center of the pharmacophoric point

          radius: Quantity:
              The new radius of the pharmacophoric point
        """
        index = self.get_picked_point_index(view)
        if index is not None and index < len(self):
            self[index].center = center
            self[index].radius = radius

    def add_point_in_picked_location(self, view, feat_name, radius):
        """ Add a pharmacophoric point in the selected atom of a ngl view.

            Parameters
              ----------
              view : nglview.NGLWidget
                  View where the pharmacophore will be added.

              feat_name: str
                  The feature type of the pharmacophoric point

              radius: Quantity:
                  The new radius of the pharmacophoric point

        """
        try:
            atom = view.picked["atom1"]
        except KeyError:
            return

        center = puw.quantity([atom["x"], atom["y"], atom["z"]], "angstroms")
        new_point = PharmacophoricPoint(feat_name, center, radius)
        self._points.append(new_point)

    def add_to_view(self, view, palette=None, opacity=0.5):
        """Add the pharmacophore representation to a view from NGLView.

           Parameters
           ----------
           view : nglview.NGLWidget
               View where the pharmacophore will be added.

           palette : dict[str, str], optional
               Dictionary with a color for each feature type.

           opacity : float
                The level of opacity of the points. Must be a number between 0 and 1.

        """
        for point in self.pharmacophoric_points:
            point.add_to_ngl_view(view, palette, opacity)

    def add_ligands_to_view(self, view):
        """ Adds the ligands to a ngl view.

            Parameters
            ----------
            view : nglview.NGLWidget
               View where the pharmacophore will be added.
        """
        for ligand in self._ligands:
            component = view.add_component(ligand)
            component.clear()
            component.add_ball_and_stick(multipleBond=True)

    def show(self, ligands=True, palette=None):
        """ Show the pharmacophore model.

        Parameters
        ----------
        ligands : bool
            Whether to show the ligands from which this pharmacophore
            was extracted from.

        palette : str or dict.
            Color palette name or dictionary. (Default: 'openpharmacophore')

        Returns
        -------
        nglview.NGLWidget
            A nglview.NGLWidget with the 'view' of the pharmacophoric model.
        """
        view = nv.NGLWidget()
        if ligands:
            self.add_ligands_to_view(view)
        self.add_to_view(view, palette=palette)
        return view

    def to_json(self, file_name):
        """ Save the pharmacophore to a json file.

            Parameters
            ----------
            file_name: str
                Name of the json file.
        """
        data = io.json_pharmacophoric_elements(self.pharmacophoric_points)
        with open(file_name, "w") as fp:
            json.dump(data, fp)

    def to_ligand_scout(self, file_name):
        """ Save the pharmacophore to a ligand scout file (pml).

            Parameters
            ----------
            file_name: str
                Name of the json file.
        """
        xml_tree = io.ligandscout_xml_tree(self.pharmacophoric_points)
        xml_tree.write(file_name, encoding="UTF-8", xml_declaration=True)

    def to_moe(self, file_name):
        """ Save the pharmacophore to a moe file (ph4).

            Parameters
            ----------
            file_name: str
                Name of the json file.
        """
        pharmacophore_str = io.ph4_string(self.pharmacophoric_points)
        with open(file_name, "w") as fp:
            fp.write(pharmacophore_str)

    def to_mol2(self, file_name):
        """ Save the pharmacophore to a mol2 file.

            Parameters
            ----------
            file_name: str
                Name of the json file.
        """
        pharmacophore_data = io.mol2_file_info([self.pharmacophoric_points])
        with open(file_name, "w") as fp:
            fp.writelines(pharmacophore_data[0])

    def to_rdkit(self):
        """ Returns a rdkit pharmacophore with the pharmacophoric_points from the original pharmacophore.

            rdkit pharmacophores do not store the pharmacophoric_points radii, so they are returned as well.

            Returns
            -------
            rdkit_pharmacophore : rdkit.Chem.Pharmacophore
                The rdkit pharmacophore.

            radii : list[float]
                List with the radius in angstroms of each pharmacophoric point.
        """
        return rdkit_pharmacophore(self.pharmacophoric_points)

    def distance_matrix(self):
        """ Compute the distance matrix of the pharmacophore.

            Returns
            -------
            dis_matrix : np.ndarray of shape(num_points, num_points)
                The distance matrix.
        """
        n_pharmacophoric_points = self.num_points
        dis_matrix = np.zeros((n_pharmacophoric_points, n_pharmacophoric_points))

        for ii in range(n_pharmacophoric_points):
            for jj in range(ii, n_pharmacophoric_points):
                if ii == jj:
                    dis_matrix[ii, jj] = 0
                else:
                    distance = distance_between_pharmacophoric_points(
                        self[ii], self[jj])
                    dis_matrix[ii, jj] = distance
                    dis_matrix[jj, ii] = distance

        return dis_matrix

    @staticmethod
    def _is_ligand_file(file_name):
        """ Check if a file belongs to a molecular file format.

            Returns
            -------
            bool
        """
        file_extension = file_name.split(".")[-1]
        if file_extension in ["mol2", "smi", "sdf"]:
            return True
        return False

    def add_hydrogens(self, ligands="all"):
        """ Add hydrogens to one or more ligands.

            Ligands must have hydrogens added prior to pharmacophore extraction.

            Parameters
            ----------
            ligands : 'all' or list[int]
                The indices of the ligands to which hydrogens will be added.
        """
        if ligands == "all":
            self._ligands = [Chem.AddHs(lig) for lig in self._ligands]
        else:
            lig_hyd = []
            for ii, lig in enumerate(self._ligands):
                if ii in ligands:
                    lig_hyd.append(Chem.AddHs(lig))
                else:
                    lig_hyd.append(lig)
            self._ligands = lig_hyd

    def generate_conformers(self, n_confs, ligands="all", random_seed=-1):
        """ Add conformers to a ligand. It is recommended to add hydrogens before
            generating conformers.

            Parameters
            ----------
            n_confs : int or list[int]
                Number of conformers to generate for each ligand. If a single number is
                given the same number of conformers will be generated for all ligands.

            ligands : 'all' or list[int]
                The indices of the ligands to which conformers will be added.

            random_seed : int, optional
                Random seed to use.
        """
        if ligands == "all":
            lig_ind = list(range(len(self._ligands)))
        else:
            lig_ind = ligands

        if isinstance(n_confs, int):
            n_confs = [n_confs] * len(lig_ind)

        if len(n_confs) != len(lig_ind):
            raise ValueError("n_confs must have the same size as ligands")

        for ii in range(len(lig_ind)):
            Chem.EmbedMultipleConfs(
                self._ligands[lig_ind[ii]], numConfs=n_confs[ii], randomSeed=random_seed)

    def extract(self, n_points, min_actives=None):
        """ Extracts and scores pharmacophores from a set of ligands.

            Parameters
            ----------
            n_points : int
                Extracted pharmacophores will have this number of pharmacophoric
                points.

            min_actives : int
                Number of ligands that must match a common pharmacophore.
        """
        if min_actives is None:
            min_actives = len(self._ligands)
        self._pharmacophores, scores = find_common_pharmacophores(self.ligands, n_points, min_actives)

    def __len__(self):
        return len(self._points)

    def __getitem__(self, index):
        return self._points[index]

    def __eq__(self, other):
        """ Check equality between pharmacophores.

            Assumes that pharmacophoric points are sorted equally in both pharmacophores.
        """
        if isinstance(other, type(self)) and self.num_points == other.num_points:
            for ii in range(self.num_points):
                if not self.pharmacophoric_points[ii] == other.pharmacophoric_points[ii]:
                    return False
            return True
        return False

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(n_pharmacophoric_points: {self.num_points})"
