from openpharmacophore import PharmacophoricPoint, Pharmacophore
from openpharmacophore.algorithms.cliques import PharmacophoreGraph
import numpy as np
import pytest
import pyunitwizard as puw


@pytest.fixture
def three_element_pharmacophore():
    """ Returns a pharmacophore with three pharmacophoric points. """
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

    return Pharmacophore(pharmacophoric_points=points)

@pytest.fixture
def four_element_pharmacophore():
    """ Returns a pharmacophore with four pharmacophoric points."""
    radius = puw.quantity(1.0, "angstroms")
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
    
    return Pharmacophore(pharmacophoric_points=points)

def test_adjacency_matrix(three_element_pharmacophore, four_element_pharmacophore):
    
    matrix_1 = PharmacophoreGraph.get_adjacency_matrix(three_element_pharmacophore, dmin=2.0, dmax=13.0, bin_size=1)
    
    assert matrix_1.shape == (3, 3)
    assert np.all(matrix_1 == np.array([["A", 5.0, 6.0],
                                        [5.0, "D", 5.0],
                                        [6.0, 5.0, "R"]],
                                        dtype=object))
    
    matrix_2 = PharmacophoreGraph.get_adjacency_matrix(four_element_pharmacophore, dmin=1.0, dmax=10.0, bin_size=1)

    assert matrix_2.shape == (4, 4)
    assert np.all(matrix_2 == np.array([["A", 4.0, 5.0, 4.0],
                                        [4.0, "N", 1.0, 2.0],
                                        [5.0, 1.0, "P", 4.0],
                                        [4.0, 2.0, 4.0, "R"]],
                                        dtype=object))
    
    matrix_2 = PharmacophoreGraph.get_adjacency_matrix(four_element_pharmacophore, dmin=1.0, dmax=10.0, bin_size=1, delta=0.5)

    assert matrix_2.shape == (4, 4)
    assert np.all(matrix_2 == np.array([["A", 4.0, 5.0, 4.0],
                                        [3.0, "N", 1.0, 2.0],
                                        [4.0, 1.0, "P", 4.0],
                                        [3.0, 2.0, 3.0, "R"]],
                                        dtype=object))