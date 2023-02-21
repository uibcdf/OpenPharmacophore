import nglview as nv
import pyunitwizard as puw
from openpharmacophore import Protein, ComplexBindingSite, Ligand, Pharmacophore
from openpharmacophore.molecular_systems.chem_feats import ChemFeatContainer
from openpharmacophore.io.pdb import write_pdb_block
from openpharmacophore import constants


class Viewer:
    """ Class to visualize molecular systems and pharmacophores.
    """

    def __init__(self):
        self._widget = nv.NGLWidget()
        self._has_protein = False
        self._has_ligand = False
        self._has_chem_feats = False
        self._has_pharmacophore = False

    @property
    def n_components(self) -> int:
        return len(self._widget._ngl_component_names)

    @property
    def has_protein(self) -> bool:
        """ Returns true if the view contains a protein."""
        return self._has_protein

    @property
    def has_ligand(self) -> bool:
        """ Returns true if the view contains a ligand. """
        return self._has_ligand

    @property
    def has_chem_feats(self) -> bool:
        """ Returns true if the view contains chemical features. """
        return self._has_chem_feats

    @property
    def has_pharmacophore(self) -> bool:
        """ Returns true if the view contains a pharmacophore. """
        return self._has_pharmacophore

    def add_components(self, components):
        """ Add multiple components to the viewer, such as Proteins, Ligands and
            BindingSites.

            Parameters
            ----------
            components : list[Any]
        """
        for comp in components:
            self.add_component(comp)

    def add_component(self, component):
        """ Add a single component to the viewer, such as a Protein, Ligand or
            BindingSite.

            Parameters
            ----------
            component : Any
        """
        if isinstance(component, Protein):
            text_struct = nv.TextStructure(write_pdb_block(component.topology, component.coords), ext="pdb")
            self._widget.add_component(text_struct)
            self._has_protein = True
        elif isinstance(component, ComplexBindingSite):
            self._widget.add_component(component.to_rdkit())
            self._has_protein = True
        elif isinstance(component, Ligand):
            self._widget.add_component(component.to_rdkit())
            self._has_ligand = True
        elif isinstance(component, ChemFeatContainer):
            self.add_chem_feats(component)
        elif isinstance(component, Pharmacophore):
            self.add_pharmacophore(component)
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
        self._widget.update_representation(
            component=self.n_components+1, repr_index=0, opacity=0.9
        )

    def _add_arrow(self, centroid, length, direction, color):
        """ Add an arrow to the view

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
        # TODO: opacity is not working
        self._widget.update_representation(component=self.n_components+1, repr_index=0, opacity=0.9)

    def show(self):
        """ Shows the view.

            Returns
            -------
            nv.NGLWidget
        """
        return self._widget

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
