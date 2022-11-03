import openpharmacophore.utils.maths as maths
import pyunitwizard as puw
import numpy as np


def test_points_distance():
    distance = maths.points_distance(
        puw.quantity(np.array([0., 0., 0.]), "angstroms"),
        puw.quantity(np.array([4., 3., 0.]), "angstroms"),
    )
    assert puw.get_value(distance) == 5.


def test_quantity_norm():
    qt = puw.quantity(np.array([4, 3, 0]), "angstroms")
    assert maths.quantity_norm(qt) == 5


def test_ring_normal():
    indices = [0, 1, 2, 3, 4, 5]
    centroid = puw.quantity(np.array([-1, 1, 2]), "angstroms")
    coords = puw.quantity(np.array([
        [-4, 2, 2],
        [-2, 1, 5],
        [2, 3, 4],
        [4, 5, 6],
        [7, 8, 9],
        [10, 11, 12]
    ]), "angstroms")
    normal = maths.ring_normal(indices, coords, centroid)
    normal_expected = np.array(
        [3 / np.sqrt(91), 9 / np.sqrt(91), 1 / np.sqrt(91)])
    assert np.allclose(normal, normal_expected)


def test_point_projection():
    normal = np.array([0, 0, 1])
    plane_point = puw.quantity(np.array([0, 0, 3]), "angstroms")
    point = puw.quantity(np.array([-5, 6, 10]), "angstroms")

    projection = maths.point_projection(normal, plane_point, point)
    proj_expected = puw.quantity((np.array([-5, 6, 3])), "angstroms")
    assert np.allclose(projection, proj_expected)


def test_angle_between_normals_vectors_equal():
    n1 = np.array([0, 0, 1])
    n2 = np.array([0, 0, 1])
    assert maths.angle_between_normals(n1, n2) == 0


def test_angle_between_normals():
    n1 = np.array([0, 0, 1])
    n2 = np.array([1, 0, 0])
    assert maths.angle_between_normals(n1, n2) == 90

    n1 = np.array([0, 0, 1])
    n2 = np.array([-1, 0, 0])
    assert maths.angle_between_normals(n1, n2) == 90


def test_maximum_distance():
    centroid = puw.quantity(np.array([0., 0., 0.]), "angstroms")
    coordinates = puw.quantity(np.array([
        [1., 1., 1.],
        [2., 2., 2.],
        [3., 3., 3.],
        [4., 4., 4.],
    ]), "angstroms")
    expected = puw.quantity(np.sqrt(48), "angstroms")
    assert maths.maximum_distance(centroid, coordinates) == expected
