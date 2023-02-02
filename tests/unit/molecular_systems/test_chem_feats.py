import openpharmacophore.molecular_systems.chem_feats as cf
import pyunitwizard as puw
import numpy as np
from rdkit import Chem


def test_create_chem_feats():
    indices = [(0, 1), (2,)]
    coords = puw.quantity(np.array([
        [1., 1., 1.],
        [2., 2., 2.],
        [3., 3., 3.],
        [4., 4., 4.],
    ]), "nanometers")
    feats = cf.create_chem_feats("hb donor", indices, coords)

    assert len(feats) == 2

    centroid_1 = puw.quantity(np.array([1.5, 1.5, 1.5]), "nanometers")
    assert feats[0].type == "hb donor"
    assert np.allclose(feats[0].coords, centroid_1)

    centroid_2 = puw.quantity(np.array([3., 3., 3.]), "nanometers")
    assert feats[1].type == "hb donor"
    assert np.allclose(feats[1].coords, centroid_2)


def test_add_features_to_container():
    feats = [
        cf.ChemFeat(type="aromatic ring", coords=puw.quantity(np.array([1.5, 1.5, 1.5]), "nanometers")),
        cf.ChemFeat(type="hb acceptor", coords=puw.quantity(np.array([2.5, 2.5, 2.5]), "nanometers")),
    ]
    container = cf.ChemFeatContainer()
    container.add_feats(feats)

    assert container.has_feat("aromatic ring")
    assert container.has_feat("hb acceptor")

    assert len(container.aromatic) == 1
    assert len(container.acceptor) == 1


def test_mol_chem_feats():
    indices = {
        "hydrophobicity": [(0, 1)],
        "positive charge": [(3, )]
    }
    coords = puw.quantity(np.array([
        [1., 1., 1.],
        [2., 2., 2.],
        [3., 3., 3.],
        [4., 4., 4.],
    ]), "nanometers")

    container = cf.mol_chem_feats(indices, coords)

    assert len(container) == 2
    assert len(container.hydrophobic) == 1
    assert len(container.positive) == 1


def test_aromatic_chem_feats_have_normal_vector():
    indices = [(0, 1, 2, 3, 4, 5)]
    coords = puw.quantity(np.array([
        # ring is an hexagon on the xy plane
        [1, 0, 0],
        [1 / 2, np.sqrt(3) / 2, 0],
        [-1 / 2, np.sqrt(3) / 2, 0],
        [-1, 0, 0],
        [-1 / 2, -np.sqrt(3) / 2, 0],
        [1 / 2, -np.sqrt(3) / 2, 0],
    ]), "angstroms")
    feats = cf.aromatic_chem_feats(indices, coords)

    assert len(feats) == 1
    normal = feats[0].normal
    centroid = feats[0].coords

    assert np.allclose(normal, np.array([0., 0., 1.]))
    assert np.allclose(centroid, np.zeros((3,)))


def test_hb_donor_chem_feats_have_hydrogen_coordinates():
    indices = [(0,)]
    water = Chem.MolFromSmiles("O")
    water = Chem.AddHs(water)
    coords = puw.quantity(np.array([
       [0., 0., 0.],
       [1., 1., 1.],
       [2., 2., 2.],
    ]), "nanometers")

    feats = cf.donor_chem_feats(indices, coords, water)
    # Oxygen donor has 2 hydrogen atoms, so we expect 2 donor with the
    # same centroid but different hydrogen atom coordinates
    assert len(feats) == 2
    assert np.all(feats[0].hyd == coords[1])
    assert np.all(feats[1].hyd == coords[2])
