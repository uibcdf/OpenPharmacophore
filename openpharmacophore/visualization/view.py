import nglview as nv
import pyunitwizard as puw
from openpharmacophore import Protein, ComplexBindingSite, Ligand, Pharmacophore
from openpharmacophore.molecular_systems.chem_feats import ChemFeatContainer
from openpharmacophore.io.pdb import write_pdb_block
from openpharmacophore import LigandBasedPharmacophore, LigandReceptorPharmacophore
from openpharmacophore import constants


class Viewer:
    """ Class to visualize molecular systems and pharmacophores.

        Wrapper around nglview NGLWidget
    """
    # TODO: method to convert to HTML

    def __init__(self):
        self._widget = nv.NGLWidget()
        self._components = []

        self._has_protein = False
        self._has_ligand = False
        self._has_chem_feats = False
        self._has_pharmacophore = False

    @property
    def n_components(self) -> int:
        return len(self._widget._ngl_component_names)

    @property
    def has_protein(self) -> bool:
        """ Returns true if the visualization contains a protein."""
        return self._has_protein

    @property
    def has_ligand(self) -> bool:
        """ Returns true if the visualization contains a ligand. """
        return self._has_ligand

    @property
    def has_chem_feats(self) -> bool:
        """ Returns true if the visualization contains chemical features. """
        return self._has_chem_feats

    @property
    def has_pharmacophore(self) -> bool:
        """ Returns true if the visualization contains a pharmacophore. """
        return self._has_pharmacophore

    def add_components(self, components):
        """ Add multiple components to the viewer, such as Proteins, Ligands and
            BindingSites.

            Parameters
            ----------
            components : list[Any]
        """
        for comp in components:
            self._components.append(comp)

    def load_component(self, component, struct):
        """ Load a component to the viewer, such as a Protein, Ligand or
            BindingSite so, it can be visualized.

            Parameters
            ----------
            component : Any

            struct : int, optional
                Structure, frame or conformer index.
        """
        if isinstance(component, Protein):
            text_struct = nv.TextStructure(
                write_pdb_block(component.topology, component.coords, conformer=0),
                ext="pdb")
            self._widget.add_component(text_struct)
            self._has_protein = True
        elif isinstance(component, ComplexBindingSite):
            self._add_molecule(component.to_rdkit(), struct)
            self._has_protein = True
        elif isinstance(component, Ligand):
            self._add_molecule(component.to_rdkit(), struct)
            self._has_ligand = True
        elif isinstance(component, ChemFeatContainer):
            self.add_chem_feats(component)
        elif isinstance(component, Pharmacophore):
            self.add_pharmacophore(component)
        elif isinstance(component, (LigandReceptorPharmacophore, LigandBasedPharmacophore)):
            self.add_pharmacophore(component[struct])
        else:
            raise NotImplementedError

    def add_chem_feats(self, chem_feats):
        """ Add chemical features to the viewer.

            Parameters
            ----------
            chem_feats : ChemFeatContainer

        """
        radius = 1.0  # in angstroms
        for feat in chem_feats:
            color = constants.PALETTE[feat.type]
            coords = puw.get_value(feat.coords, "angstroms").tolist()
            self._add_sphere(coords, radius, color, feat.type)
        self._has_chem_feats = True

    def add_pharmacophore(self, pharmacophore):
        """ Add a pharmacophore to the viewer

            Parameters
            ----------
            pharmacophore : Pharmacophore

        """
        for point in pharmacophore:
            coords = puw.get_value(point.center, "angstroms")
            color = constants.PALETTE[point.feature_name]
            radius = puw.get_value(point.radius, "angstroms")
            self._add_sphere(coords.tolist(), radius, color, point.feature_name)
            if point.has_direction:
                self._add_arrow(coords, radius, point.direction, color)
        self._has_pharmacophore = True

    def _add_sphere(self, centroid, radius, color, label):
        """ Add a sphere to the viewer.

            Parameters
            ----------
            centroid : list[float]
                A list of length 3 with the coordinates of the centroid

            radius : float
                Radius of the sphere

            color : tuple[float]
                Color in rgb format.

            label : str
                The label of the sphere
        """
        self._widget.shape.add_sphere(centroid, color, radius, label)

    def _add_arrow(self, centroid, length, direction, color):
        """ Add an arrow to the visualization

            Parameters
            ----------
            centroid : np.ndarray
                Aan array of length 3 with the coordinates of the centroid.

            length : float
                Length of the arrow.

            direction : np.ndarray
                Directionality of the arrow

            color : tuple[float]
                Color in rgb format.

        """
        arrow_radius = 0.2
        end_arrow = (centroid + length * direction * 2.0).tolist()

        self._widget.shape.add_arrow(centroid, end_arrow, color, arrow_radius)

    def _add_molecule(self, mol, conformer):
        """ Add an rdkit molecule to the viewer.

            Parameters
            ----------
            mol : rdkit.Chem.Mol

            conformer : int
        """
        self._widget.add_component(
            nv.RdkitStructure(mol, conf_id=conformer)
        )

    def _restore_widget(self):
        """ Restore the widget to have a clean visualization."""
        if self.n_components > 0:
            self._widget = nv.NGLWidget()

    def show(self, struct=0):
        """ Shows the visualization.

            Parameters
            ----------
            struct : int, optional
                Frame, structure or conformer index to show if there are any proteins,
                ligands or pharmacophores with multiple structures.

            Returns
            -------
            nv.NGLWidget
        """
        self._restore_widget()
        for comp in self._components:
            self.load_component(comp, struct)
        self._update_opacity()
        return self._widget

    def _update_opacity(self):
        """ Add opacity to chemical features and pharmacophoric points. """
        for ii, comp in enumerate(self._widget._ngl_component_names):
            if comp == "nglview.shape.Shape":
                self._widget.update_representation(component=ii, repr_index=0, opacity=0.6)

    def set_protein_style(self, style):
        """ Set the style of the protein.

            Parameters
            ----------
            style : str
                Name of the style
        """
        styles = [
            "cartoon",
            "licorice",
            "ball+stick",
        ]
        if style not in styles:
            raise ValueError
        self._widget.clear_representations()
        self._widget.add_representation(style, selection="protein")

    def to_nglview(self):
        """ Get an NGLView widget representing the view.

            Returns
            -------
            nv.NGLWidget
        """
        return self._widget
