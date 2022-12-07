from openpharmacophore import PharmacophoricPoint, Pharmacophore
from openpharmacophore.pharmacophore.ligand_based.rdp import find_common_pharmacophores
from openpharmacophore.pharmacophore.rdkit_pharmacophore import rdkit_pharmacophore
from openpharmacophore.pharmacophore.chem_feats import smarts_ligand, feature_indices
from openpharmacophore.utils.conformers import ConformerGenerator
import openpharmacophore.io as io
from openpharmacophore._private_tools.exceptions import InvalidFileFormat
from matplotlib.colors import to_rgb
import nglview as nv
import pyunitwizard as puw
import rdkit.Chem.AllChem as Chem
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG
from copy import deepcopy
import json
import os
import math


class LigandBasedPharmacophore:
    """ Class to store and extract pharmacophores from a set of ligands.

        Attributes
        ----------
        pharmacophores : list[Pharmacophore]
            List with common pharmacophores.
    """

    def __init__(self):
        self._pharmacophores = []  # type: list[Pharmacophore]
        self._ligands = []
        self._feats = []

    @property
    def ligands(self):
        return self._ligands

    @ligands.setter
    def ligands(self, ligand_list):
        self._ligands = ligand_list

    @ligands.deleter
    def ligands(self):
        self._ligands.clear()

    @property
    def pharmacophores(self):
        return self._pharmacophores

    @property
    def feats(self):
        return self._feats

    def add_pharmacophore(self, pharma):
        """ Add a new pharmacophore.

            Parameters
            ----------
            pharma : Pharmacophore
        """
        self._pharmacophores.append(pharma)

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
        self._pharmacophores.append([])
        if file_name.endswith(".json"):
            self._pharmacophores[0] = io.load_json_pharmacophore(file_name)[0]
        elif file_name.endswith(".mol2"):
            self._pharmacophores[0] = io.load_mol2_pharmacophoric_points(file_name)[0]
        elif file_name.endswith(".pml"):
            self._pharmacophores[0] = io.read_ligandscout(file_name)
        elif file_name.endswith(".ph4"):
            self._pharmacophores[0] = io.pharmacophoric_points_from_ph4_file(file_name)
        else:
            raise InvalidFileFormat(file_name.split(".")[-1])

    def add_point(self, point, pharma):
        """ Adds a pharmacophoric point.

            Parameters
            ----------
            point : PharmacophoricPoint
                The pharmacophoric point that will be added.

            pharma : int
                Index of the pharmacophore to which the point will be added
        """
        self._pharmacophores[pharma].add(point)

    def remove_point(self, pharma, point):
        """ Removes a pharmacophoric point from the pharmacophore.

            Parameters
            ----------
            pharma : int
                Index of the pharmacophore.

            point : int
                The index of the pharmacophoric point.
        """
        self._pharmacophores[pharma].remove(point)

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

    def remove_picked_point(self, view, pharma):
        """ Remove a pharmacophoric point selected in a ngl view.

            If nothing is selected or the selected thing is not a
            pharmacophoric point, it doesn't change the view.

           Parameters
           ----------
           view : nglview.NGLWidget
               View where the pharmacophore will be added.

            pharma : int
                Index of pharmacophore
        """
        index = self.get_picked_point_index(view)
        if index is not None and index < len(self):
            self.remove_point(pharma, index)

    def edit_picked_point(self, view, center, radius, pharma):
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

           pharma : int
                Index of pharmacophore
        """
        index = self.get_picked_point_index(view)
        if index is not None and index < len(self):
            self[pharma][index].center = center
            self[pharma][index].radius = radius

    def add_point_in_picked_location(self, view, feat_name, radius, pharma):
        """ Add a pharmacophoric point in the selected atom of a ngl view.

            Parameters
              ----------
              view : nglview.NGLWidget
                  View where the pharmacophore will be added.

              feat_name: str
                  The feature type of the pharmacophoric point

              radius: Quantity:
                  The new radius of the pharmacophoric point.

              pharma : int
                Index of pharmacophore

        """
        try:
            atom = view.picked["atom1"]
        except KeyError:
            return

        center = puw.quantity([atom["x"], atom["y"], atom["z"]], "angstroms")
        new_point = PharmacophoricPoint(feat_name, center, radius)
        self._pharmacophores[pharma].add(new_point)

    def add_to_view(self, view, pharma=0, palette=None, opacity=0.5):
        """Add the pharmacophore representation to a view from NGLView.

           Parameters
           ----------
           view : nglview.NGLWidget
               View where the pharmacophore will be added.

            pharma : int
                Index of the pharmacophore

           palette : dict[str, str], optional
               Dictionary with a color for each feature type.

           opacity : float
                The level of opacity of the points. Must be a number between 0 and 1.

        """
        for point in self[pharma]:
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
        # TODO: save multiple pharmacophores
        data = io.json_pharmacophoric_elements(self.pharmacophores[0])
        with open(file_name, "w") as fp:
            json.dump(data, fp)

    def to_ligand_scout(self, file_name):
        """ Save the pharmacophore to a ligand scout file (pml).

            Parameters
            ----------
            file_name: str
                Name of the json file.
        """
        xml_tree = io.ligandscout_xml_tree(self.pharmacophores[0])
        xml_tree.write(file_name, encoding="UTF-8", xml_declaration=True)

    def to_moe(self, file_name):
        """ Save the pharmacophore to a moe file (ph4).

            Parameters
            ----------
            file_name: str
                Name of the json file.
        """
        pharmacophore_str = io.ph4_string(self.pharmacophores[0])
        with open(file_name, "w") as fp:
            fp.write(pharmacophore_str)

    def to_mol2(self, file_name):
        """ Save the pharmacophore to a mol2 file.

            Parameters
            ----------
            file_name: str
                Name of the json file.
        """
        pharmacophore_data = io.mol2_file_info([self.pharmacophores[0]])
        with open(file_name, "w") as fp:
            fp.writelines(pharmacophore_data[0])

    def to_rdkit(self, pharma):
        """ Returns a rdkit pharmacophore with the pharmacophoric_points from the original pharmacophore.

            Returns
            -------
            rdkit_pharmacophore : rdkit.Chem.Pharmacophore
                The rdkit pharmacophore.

            pharma : int
                Index of the pharmacophore.
        """
        return rdkit_pharmacophore(self.pharmacophores[pharma])

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

    def generate_conformers(self, n_confs=-1, ligands="all",  *args, **kwargs):
        """ Add conformers to ligands. It is recommended to add hydrogens before
            generating conformers.

            Parameters
            ----------
            n_confs : int or list[int]
                Number of conformers to generate for each ligand. If a single number is
                given the same number of conformers will be generated for all ligands.

            ligands : 'all' or list[int]
                The indices of the ligands to which conformers will be added.

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
            conf_gen = ConformerGenerator(n_confs[ii], *args, **kwargs)
            self._ligands[lig_ind[ii]] = conf_gen(self._ligands[lig_ind[ii]])

    def find_chem_feats(self):
        """ Find chemical features in the ligands associated with this pharmacophore.

            The attribute feats is update with a dictionary containing the chemical
            features of each ligand.
        """
        for ii, lig in enumerate(self._ligands):
            self._feats.append({})
            for feat_type, smarts in smarts_ligand.items():
                indices = feature_indices(smarts, lig)
                short_name = PharmacophoricPoint.feature_to_char[feat_type]
                self._feats[ii][short_name] = indices

    def extract(self, n_points, min_actives=None, max_pharmacophores=None):
        """ Extracts and scores pharmacophores from a set of ligands.

            Parameters
            ----------
            n_points : int
                Extracted pharmacophores will have this number of pharmacophoric
                points.

            min_actives : int, optional
                Number of ligands that must match a common pharmacophore.

            max_pharmacophores : int, optional
                Maximum number of pharmacophores to return. If set to null
                all found pharmacophores will be returned.
        """
        if min_actives is None:
            min_actives = len(self._ligands)
        self._pharmacophores = find_common_pharmacophores(
            self.ligands, self.feats, n_points, min_actives, max_pharmacophores
        )

    def draw(self, mol_size):
        """ Draw the ligands with their chemical features highlighted.

            Parameters
            ---------
            mol_size : tuple[int, int]
                A tuple with the width and height of each ligand drawing.

            Returns
            -------
            drawing
        """
        ligands = [deepcopy(lig) for lig in self.ligands]
        for lig in ligands:
            lig.RemoveAllConformers()

        width, height = self._drawing_size(mol_size[0], mol_size[1], len(self._ligands))
        atoms, colors, radii = self._atom_highlights()
        drawing = MolDraw2DSVG(width=width, height=height,
                               panelWidth=mol_size[0], panelHeight=mol_size[1])
        drawing.DrawMolecules(ligands,
                              highlightAtoms=atoms,
                              highlightAtomColors=colors,
                              highlightAtomRadii=radii
                              )
        drawing.FinishDrawing()
        return drawing

    @staticmethod
    def _drawing_size(mol_width, mol_height, n_mols):
        """ Returns the width and height of a drawing consisting of n_mols
            with the given width and height

            Parameters
            ----------
            mol_width : int
            mol_height : int
            n_mols : int

            Returns
            -------
            width : int
            height : int
        """
        max_per_row = 4
        max_width = max_per_row * mol_width
        if mol_width * n_mols <= max_width:
            width = mol_width * n_mols
            n_rows = 1
        else:
            width = max_width
            n_rows = math.ceil(mol_width * n_mols / width)
        height = n_rows * mol_height
        return width, height

    def _atom_highlights(self):
        """ Returns the indices, colors and radii of the atoms that will be
            highlighted in a drawing

            Returns
            -------
            atoms : list[tuple[int]]
            colors : list[dict[int, tuple[int]]
            radii : list[dict[int, float]]
        """
        atoms = []
        colors = []
        radii = []

        for lig_feats in self._feats:
            lig_ind = []
            lig_color = {}
            lig_radii = {}
            for feat_type, ind_list in lig_feats.items():
                color = PharmacophoricPoint.palette[PharmacophoricPoint.char_to_feature[feat_type]]
                color = to_rgb(color)
                for ind_tup in ind_list:
                    for ind in ind_tup:
                        lig_ind.append(ind)
                        lig_color[ind] = color
                        lig_radii[ind] = 0.5
            atoms.append(tuple(lig_ind))
            colors.append(lig_color)
            radii.append(lig_radii)

        return atoms, colors, radii

    def __len__(self):
        return len(self._pharmacophores)

    def __getitem__(self, index):
        return self._pharmacophores[index]

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(n_pharmacophores: {len(self)})"
