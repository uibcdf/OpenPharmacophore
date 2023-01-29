from openpharmacophore import Viewer
from openpharmacophore.molecular_systems.chem_feats import ChemFeatContainer, ChemFeat
import pyunitwizard as puw
import numpy as np


class FakeShape:

    def add_sphere(self, centroid, color, radius, label):
        pass


class FakeWidget:

    def __init__(self):
        self._ngl_component_names = []
        self.shape = FakeShape()

    def add_component(self, component):
        self._ngl_component_names.append(component)

    def update_representation(self, component, repr_index, opacity):
        self._ngl_component_names.append(component)


def get_viewer():
    viewer = Viewer()
    viewer._widget = FakeWidget()
    return viewer


def test_add_components_to_viewer(
        protein_4_residues, estradiol
):
    viewer = get_viewer()
    viewer.add_components(
        [protein_4_residues, estradiol]
    )
    assert viewer.n_components == 2
    assert viewer.has_protein
    assert viewer.has_ligand


def test_add_chem_feats_to_viewer():
    viewer = get_viewer()
    chem_feats = ChemFeatContainer(
        [ChemFeat("aromatic ring", puw.quantity(np.array([1.5, 1.5, 1.5]), "nanometers")),
         ChemFeat("hb acceptor", puw.quantity(np.array([2.5, 2.5, 2.5]), "nanometers"))]
    )
    viewer.add_chem_feats(chem_feats)

    assert viewer.n_components == 2
    assert viewer.has_chem_feats
