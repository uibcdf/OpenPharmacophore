from openpharmacophore import Pharmacophore, PharmacophoricPoint
from openpharmacophore._private_tools.exceptions import InvalidFeatureError
import pyunitwizard as puw
import pytest
import networkx as nx
import numpy as np
from rdkit.Chem import Pharm3D
import copy

@pytest.fixture
def three_point_pharmacophore() -> Pharmacophore:
    """Returns as pharmacophore with three pharmacophoric_points"""
    hb_acceptor = PharmacophoricPoint(
        feat_type="hb acceptor",
        center=puw.quantity([1,0,0], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    ring_1 = PharmacophoricPoint(
        feat_type="aromatic ring",
        center=puw.quantity([2, 1, 4], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    ring_2 = PharmacophoricPoint(
        feat_type="aromatic ring",
        center=puw.quantity([0, 1, 2], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    return Pharmacophore(pharmacophoric_points=[ring_1, hb_acceptor, ring_2])

@pytest.fixture
def two_point_pharmacophore() -> Pharmacophore:
    """Returns a pharmacophore with two pharmacophoric points."""
    hb_donor = PharmacophoricPoint(
        feat_type="hb donor",
        center=puw.quantity([1,0,1], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    hydrophobicity = PharmacophoricPoint(
        feat_type="hydrophobicity",
        center=puw.quantity([-1, 0, 3], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    return Pharmacophore(pharmacophoric_points=[hb_donor, hydrophobicity])

@pytest.fixture
def five_point_pharmacophore(three_point_pharmacophore) -> Pharmacophore:
    hydrophobicity = PharmacophoricPoint(
        feat_type="hydrophobicity",
        center=puw.quantity([-1, 0, 3], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    positive =  PharmacophoricPoint(
        feat_type="positive charge",
        center=puw.quantity([3, -1, 4], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    pharmacophore = copy.deepcopy(three_point_pharmacophore)
    pharmacophore.add_point(hydrophobicity)
    pharmacophore.add_point(positive)
    return pharmacophore

def test_get_point(three_point_pharmacophore):
    # Test that we can get an element by indexing and calling the get_point method
    assert three_point_pharmacophore[0] == three_point_pharmacophore.get_point(0)
    assert three_point_pharmacophore[1] == three_point_pharmacophore.get_point(1)
    assert three_point_pharmacophore[2] == three_point_pharmacophore.get_point(2)
    assert three_point_pharmacophore[0].feature_name == "hb acceptor"
    assert three_point_pharmacophore[1].feature_name == "aromatic ring"
    assert three_point_pharmacophore[2].feature_name == "aromatic ring"
    
def test_inmutable(three_point_pharmacophore):
    # Test that we cannot modify the pharmacophoric point list
    hydrophobicity = PharmacophoricPoint(
        feat_type="hydrophobicity",
        center=puw.quantity([-1, 0, 3], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )

    with pytest.raises(TypeError):
        three_point_pharmacophore[2] = hydrophobicity

def test_iter(three_point_pharmacophore):
    for point in three_point_pharmacophore:
        assert isinstance(point, PharmacophoricPoint)

def test_equality(three_point_pharmacophore, two_point_pharmacophore):
    
    other_pharmacophore = copy.deepcopy(three_point_pharmacophore)
    assert three_point_pharmacophore == other_pharmacophore
    assert three_point_pharmacophore != two_point_pharmacophore

def test_representation(three_point_pharmacophore, two_point_pharmacophore):
    assert three_point_pharmacophore.__repr__() == "Pharmacophore(n_pharmacophoric_points: 3)"
    assert two_point_pharmacophore.__repr__() == "Pharmacophore(n_pharmacophoric_points: 2)"

def test_add_point(three_point_pharmacophore):

    hydrophobic = PharmacophoricPoint(
        feat_type="hydrophobicity",
        center=puw.quantity([-1, 0, 2], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )

    positive =  PharmacophoricPoint(
        feat_type="positive charge",
        center=puw.quantity([3, -1, 4], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )

    pharmacophore = copy.deepcopy(three_point_pharmacophore)
    pharmacophore.add_point(hydrophobic)
    assert pharmacophore.n_pharmacophoric_points == len(pharmacophore)
    assert pharmacophore.n_pharmacophoric_points == 4
    # Element should be sorted
    assert pharmacophore[1].feature_name == "hydrophobicity"

    pharmacophore.add_point(positive)
    assert pharmacophore.n_pharmacophoric_points == len(pharmacophore)
    assert pharmacophore.n_pharmacophoric_points == 5
    assert pharmacophore[2].feature_name == "positive charge"

    # Test adding a repetead feature type
    hb_acceptor = PharmacophoricPoint(
        feat_type="hb acceptor",
        center=puw.quantity([-1,0,-1], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    pharmacophore.add_point(hb_acceptor)
    assert pharmacophore.n_pharmacophoric_points == len(pharmacophore)
    assert pharmacophore.n_pharmacophoric_points == 6
    assert pharmacophore[0].feature_name == "hb acceptor"
    assert pharmacophore[1].feature_name == "hb acceptor"

    ring = PharmacophoricPoint(
        feat_type="aromatic ring",
        center=puw.quantity([-5, -1, -3], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )

    pharmacophore.add_point(ring)
    assert pharmacophore.n_pharmacophoric_points == len(pharmacophore)
    assert pharmacophore.n_pharmacophoric_points == 7
    assert pharmacophore[4].feature_name == "aromatic ring"
    assert pharmacophore[5].feature_name == "aromatic ring"
    assert pharmacophore[6].feature_name == "aromatic ring"

def test_remove_points(three_point_pharmacophore, five_point_pharmacophore):

    # Test pharmacophore with three elements
    pharmacophore = copy.deepcopy(three_point_pharmacophore)
   

    with pytest.raises(IndexError):
        pharmacophore.remove_points([3, 4])

    pharmacophore.remove_point(0)
    assert len(pharmacophore) == 2
    assert pharmacophore.n_pharmacophoric_points == 2
    assert pharmacophore[0].feature_name == "aromatic ring"
    assert pharmacophore[1].feature_name == "aromatic ring"

    pharmacophore.remove_points([0,1])
    assert len(pharmacophore) == 0
    assert pharmacophore.n_pharmacophoric_points == 0

    # Test pharmacophore with five elements
    pharmacophore = copy.deepcopy(five_point_pharmacophore)
    pharmacophore.remove_points([1, 4])
    assert len(pharmacophore) == 3
    assert pharmacophore.n_pharmacophoric_points == 3
    assert pharmacophore[0].feature_name == "hb acceptor"
    assert pharmacophore[1].feature_name == "positive charge"
    assert pharmacophore[2].feature_name == "aromatic ring"

    pharmacophore.remove_point(1)
    assert len(pharmacophore) == 2
    assert pharmacophore.n_pharmacophoric_points == 2
    assert pharmacophore[0].feature_name == "hb acceptor"
    assert pharmacophore[1].feature_name == "aromatic ring"

    pharmacophore.remove_point(0)
    assert len(pharmacophore) == 1
    assert pharmacophore.n_pharmacophoric_points == 1
    assert pharmacophore[0].feature_name == "aromatic ring"

    pharmacophore.remove_point(0)
    assert len(pharmacophore) == 0
    assert pharmacophore.n_pharmacophoric_points == 0


def test_remove_feature(three_point_pharmacophore, five_point_pharmacophore):

    pharmacophore = copy.deepcopy(three_point_pharmacophore)
    pharmacophore.remove_feature("aromatic ring")
    assert len(pharmacophore) == 1
    assert pharmacophore[0].feature_name == "hb acceptor"

    pharmacophore.remove_feature("hb acceptor")
    assert len(pharmacophore) == 0

    pharmacophore = copy.deepcopy(five_point_pharmacophore)
    pharmacophore.remove_feature("aromatic ring")
    assert len(pharmacophore) == 3
    assert pharmacophore[0].feature_name == "hb acceptor"
    assert pharmacophore[1].feature_name == "hydrophobicity"
    assert pharmacophore[2].feature_name == "positive charge"

def test_reset(three_point_pharmacophore):
    pharmacophore = copy.deepcopy(three_point_pharmacophore)
    pharmacophore._reset()
    assert pharmacophore.n_pharmacophoric_points == 0
    assert len(pharmacophore._pharmacophoric_points) == 0

def test_feature_count(three_point_pharmacophore, five_point_pharmacophore):
    
    counter = {
            "aromatic ring":   2,
            "hydrophobicity":  0,
            "hb acceptor":     1,
            "hb donor":        0,
            "positive charge": 0,
            "negative charge": 0,
    }

    assert three_point_pharmacophore.feature_count() == counter

    counter = {
            "aromatic ring":   2,
            "hydrophobicity":  1,
            "hb acceptor":     1,
            "hb donor":        0,
            "positive charge": 1,
            "negative charge": 0,
    }

    assert five_point_pharmacophore.feature_count() == counter

def test_to_rdkit(three_point_pharmacophore):
    rdkit_ph, radii = three_point_pharmacophore.to_rdkit()

    assert len(radii) == 3
    assert np.allclose(np.array(radii), np.array([1.0, 1.0, 1.0]))

    assert isinstance(rdkit_ph, Pharm3D.Pharmacophore.Pharmacophore)
    feats = rdkit_ph.getFeatures()
    assert len(feats) == 3

    acceptor = feats[0]
    assert acceptor.GetFamily() == "Acceptor"
    assert np.allclose(acceptor.GetPos().x, 1.0)
    assert np.allclose(acceptor.GetPos().y, 0.0)
    assert np.allclose(acceptor.GetPos().z, 0.0)
    
    ring_1 = feats[1]
    assert ring_1.GetFamily() == "Aromatic"
    assert np.allclose(ring_1.GetPos().x, 2.0)
    assert np.allclose(ring_1.GetPos().y, 1.0)
    assert np.allclose(ring_1.GetPos().z, 4.0)

    ring_2 = feats[2]
    assert ring_2.GetFamily() == "Aromatic"
    assert np.allclose(ring_2.GetPos().x, 0.0)
    assert np.allclose(ring_2.GetPos().y, 1.0)
    assert np.allclose(ring_2.GetPos().z, 2.0)

def test_distance_matrix():

    # Test pharmacophore with zero elements
    pharmacophore = Pharmacophore()
    matrix = pharmacophore.distance_matrix()
    assert matrix.shape == (0, 0)

    radius = puw.quantity(1.0, "angstroms")
    points = [
        PharmacophoricPoint(feat_type="hb acceptor",
                            center=puw.quantity([1, 0, -5], "angstroms"),
                            radius=radius),
        PharmacophoricPoint(feat_type="hb donor",
                            center=puw.quantity([2, 1, 0], "angstroms"),
                            radius=radius),
        PharmacophoricPoint(feat_type="aromatic ring",
                            center=puw.quantity([-3, 2, -1], "angstroms"),
                            radius=radius),
    ]
    
    pharmacophore = Pharmacophore(pharmacophoric_points=points)
    distance_matrix = pharmacophore.distance_matrix()
    
    assert distance_matrix.shape == (3, 3)
    assert np.allclose(distance_matrix, 
                       np.array([[0, np.sqrt(27), 6],
                                 [np.sqrt(27), 0, np.sqrt(27)],
                                 [6, np.sqrt(27), 0]])
                       )
    
    points = [
        PharmacophoricPoint(feat_type="hb acceptor",
                            center=puw.quantity([1, 2, 3], "angstroms"),
                            radius=radius),
        PharmacophoricPoint(feat_type="negative charge",
                            center=puw.quantity([0, 0, 0], "angstroms"),
                            radius=radius),
        PharmacophoricPoint(feat_type="positive charge",
                            center=puw.quantity([-1, 0, -1], "angstroms"),
                            radius=radius),
        PharmacophoricPoint(feat_type="aromatic ring",
                            center=puw.quantity([2, -1, 1], "angstroms"),
                            radius=radius),
    ]
    
    sq = np.sqrt
    pharmacophore = Pharmacophore(pharmacophoric_points=points)
    distance_matrix = pharmacophore.distance_matrix()
    assert distance_matrix.shape == (4, 4)
    assert np.allclose(distance_matrix,
                       np.array([[0, sq(14), sq(24), sq(14)],
                                 [sq(14), 0, sq(2), sq(6)],
                                 [sq(24), sq(2), 0, sq(14)],
                                 [sq(14), sq(6), sq(14), 0]]))

def test_to_nx_graph(three_point_pharmacophore, two_point_pharmacophore, five_point_pharmacophore):
    
    # Test pharmacophore with zero elements
    pharmacophore = Pharmacophore()
    graph = pharmacophore.to_nx_graph()
    assert graph.number_of_nodes() == 0

    graph = two_point_pharmacophore.to_nx_graph()
    assert graph.number_of_nodes() == 2
    assert graph.number_of_edges() == 1

    graph = three_point_pharmacophore.to_nx_graph()
    assert isinstance(graph, nx.Graph)
    assert graph.number_of_nodes() == 3
    assert graph.number_of_edges() == 3

    edges = list(graph.edges)
    assert ("A1", "R1") in edges
    assert ("A1", "R2") in edges
    assert ("R1", "R2") in edges
    
    assert graph["A1"]["R1"]["dis"] == 4.0
    assert graph["A1"]["R2"]["dis"] == 2.0
    assert graph["R1"]["R2"]["dis"] == 3.0
    
    graph = five_point_pharmacophore.to_nx_graph()
    assert graph.number_of_nodes() == 5
    assert graph.number_of_edges() == 10