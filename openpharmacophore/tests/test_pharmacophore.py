from openpharmacophore import Pharmacophore, PharmacophoricPoint
from openpharmacophore._private_tools.exceptions import InvalidFeatureError
import pyunitwizard as puw
import pytest
import networkx as nx
import numpy as np
from rdkit.Chem import Pharm3D


@pytest.fixture
def three_point_pharmacophore():
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
    pharmacophore = Pharmacophore(pharmacophoric_points=[ring_1, hb_acceptor, ring_2])
    return pharmacophore

def test_get_point(three_element_phamracophore):
    # Test that we can get an element by indexing and calling the get_point method

    # Test that wa cannot modify the pharmacophoric point list
    pass

def test_iter(three_element_phamracophore):
    pass

def test_equality(three_element_phamrmacophore):
    pass

def test_representation(three_elment_phamrmacophore):
    pass

def test_add_element(three_point_pharmacophore):

    hydrophobic = PharmacophoricPoint(
        feat_type="hydrophobicity",
        center=puw.quantity([-1, 0, 2], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    three_point_pharmacophore.add_element(hydrophobic)
    assert three_point_pharmacophore.n_pharmacophoric_points == 4
    # Element should be sorted
    assert three_point_pharmacophore.pharmacophoric_points[1].feature_name == "hydrophobicity"

def test_remove_points(three_point_pharmacophore):
    three_point_pharmacophore.remove_pharmacophoric_points(0)
    assert three_point_pharmacophore.n_pharmacophoric_points == 2

@pytest.mark.parametrize('feat_type,exception', [
    ("aromatic ring", None), ("hydrophobicity", InvalidFeatureError)
])
def test_remove_feature(three_point_pharmacophore, feat_type, exception):
    try:
        three_point_pharmacophore.remove_feature(feat_type)
    except(InvalidFeatureError) as e:
        assert e.message == "Cannot remove feature. The pharmacophore does not contain any hydrophobicity"
    else:
        assert three_point_pharmacophore.n_pharmacophoric_points == 1
        assert three_point_pharmacophore.pharmacophoric_points[0].feature_name != "aromatic ring"

def test_reset(three_point_pharmacophore):
    three_point_pharmacophore._reset()
    assert three_point_pharmacophore.n_pharmacophoric_points == 0
    assert len(three_point_pharmacophore.pharmacophoric_points) == 0
    assert three_point_pharmacophore.extractor is None
    assert three_point_pharmacophore.molecular_system is None

def test_feature_count():
    pass

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

def test_to_nx_graph(three_point_pharmacophore):
    
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
    