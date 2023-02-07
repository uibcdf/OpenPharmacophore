import nglview as nv
import pyunitwizard as puw
from openpharmacophore import Protein, ComplexBindingSite, Ligand
from openpharmacophore.molecular_systems.convert import create_traj
from openpharmacophore.molecular_systems.chem_feats import ChemFeatContainer
from openpharmacophore import constants


class Viewer:

    def __init__(self):
        self._widget = nv.NGLWidget()
        self._has_protein = False
        self._has_ligand = False
        self._has_chem_feats = False

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
            # TODO: add Protein component without converting it to
            #   mdtraj.Trajectory
            self._widget.add_component(create_traj(
                component.coords, component.topology
            ))
            self._has_protein = True
        elif isinstance(component, ComplexBindingSite):
            self._widget.add_component(component.to_rdkit())
            self._has_protein = True
        elif isinstance(component, Ligand):
            self._widget.add_component(component.to_rdkit())
            self._has_ligand = True
        elif isinstance(component, ChemFeatContainer):
            self.add_chem_feats(component)
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
