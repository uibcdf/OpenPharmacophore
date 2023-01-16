from openpharmacophore import ComplexBindingSite
import pyunitwizard as puw
import numpy as np
import pytest


@pytest.fixture
def ligand_coords():
    return puw.quantity(
        np.array([
            [1., 2., 3.],
            [4., 5., 6.],
            [7., 8., 9.],
        ]),
        "nanometers"
    )


def test_ligand_centroid(ligand_coords):
    centroid = ComplexBindingSite._ligand_centroid(ligand_coords)
    expected = puw.quantity([4., 5., 6.], "nanometers")
    assert np.all(centroid == expected)


def test_ligand_max_extent(ligand_coords):
    centroid = puw.quantity([4., 5., 6.], "nanometers")
    max_extent = ComplexBindingSite._ligand_max_extent(ligand_coords, centroid)
    expected = puw.quantity(np.sqrt(3) * 3, "nanometers")
    assert np.allclose(max_extent, expected)
