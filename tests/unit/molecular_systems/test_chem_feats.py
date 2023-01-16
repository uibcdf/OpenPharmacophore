from openpharmacophore.molecular_systems.chem_feats import ChemFeat, ChemFeatContainer, create_chem_feats
import pyunitwizard as puw
import numpy as np


def test_create_chem_feats():
    indices = [(0, 1), (2,)]
    coords = puw.quantity(np.array([
        [1., 1., 1.],
        [2., 2., 2.],
        [3., 3., 3.],
        [4., 4., 4.],
    ]), "nanometers")
    feats = create_chem_feats("hb donor", indices, coords)

    assert len(feats) == 2

    centroid_1 = puw.quantity(np.array([1.5, 1.5, 1.5]), "nanometers")
    assert feats[0].type == "hb donor"
    assert np.allclose(feats[0].coords, centroid_1)

    centroid_2 = puw.quantity(np.array([3., 3., 3.]), "nanometers")
    assert feats[1].type == "hb donor"
    assert np.allclose(feats[1].coords, centroid_2)


def test_add_features_to_container():
    feats = [
        ChemFeat("aromatic ring", puw.quantity(np.array([1.5, 1.5, 1.5]), "nanometers")),
        ChemFeat("hb acceptor", puw.quantity(np.array([2.5, 2.5, 2.5]), "nanometers")),
    ]
    container = ChemFeatContainer()
    container.add_feats(feats)

    assert container.has_feat("aromatic ring")
    assert container.has_feat("hb acceptor")

    assert len(container.aromatic) == 1
    assert len(container.acceptor) == 1
