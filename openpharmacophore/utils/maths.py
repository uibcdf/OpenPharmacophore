import numpy as np
import pyunitwizard as puw


def points_distance(coords_1, coords_2):
    """ Returns the distance between two points in 3D space.

        Parameters
        ----------
        coords_1 : puw.Quantity
        coords_2 : puw.Quantity

        Returns
        -------
        puw.Quantity
    """
    return np.sqrt(np.sum(np.power(coords_1 - coords_2, 2)))


def quantity_norm(quantity):
    """ Returns the norm of a quantity vector.

        Parameters
        ----------
        quantity : puw.Quantity

        Returns
        -------
        float
    """
    return np.linalg.norm(puw.get_value(quantity))


def ring_normal(indices, coords, centroid):
    """ Calculate the normal vector of the plane defined by an
        aromatic ring.

        Parameters
        ----------
        indices : list[int]
            Atom indices of the ring in the trajectory.

        coords : puw.Quantity
            Coordinates array of shape=(n_atoms, 3)

        centroid : puw.Quantity
            Centroid of the ring. Array of shape=(3,)

        Returns
        -------
        normal : np.ndarray
            Vector normal to the plane. Array of shape (3,)
    """
    # Define two vectors in the plane with origin in the centroid
    vec_1 = coords[indices[0]] - centroid
    vec_2 = coords[indices[1]] - centroid
    # Normal does not have units
    normal = np.cross(puw.get_value(vec_1), puw.get_value(vec_2))
    return normal / np.linalg.norm(normal)


def point_projection(normal, plane_point, point):
    """ Returns the projection of a point into a plane.

        Parameters
        ----------
        normal : puw.Quantity
            Vector normal to the plane. Shape (3,)
        plane_point : puw.Quantity
            A point of the plane. Shape (3,)
        point : puw.Quantity
            The point that will be projected. Shape (3,)

        Returns
        -------
        puw.Quantity
            The projection of the point. Shape (3,)
    """
    # Vector from the plane to the point of interest
    vec = point - plane_point
    dist = np.dot(vec, normal)
    return point - dist * normal


def angle_between_normals(normal_1, normal_2):
    """ Compute the angle between the normals of two planes.

        Parameters
        ----------
        normal_1 : np.array
            Vector of shape (3,)

        normal_2 : np.array
            Vector of shape (3,)

        Returns
        -------
        angle: float
            Angle in degrees
    """
    denominator = np.linalg.norm(normal_1) * np.linalg.norm(normal_2)
    angle = np.arccos(np.dot(normal_1, normal_2) / denominator)
    angle = np.degrees(angle)
    if not 180 - angle < 0:
        return min(angle, 180 - angle)
    return angle


def maximum_distance(centroid, coords):
    """ Get the maximum distance from the centroid to the given
        coordinates

        Parameters
        ----------
        centroid : puw.Quantity
            Shape (3, )
        coords : puw.Quantity
            Shape (n_atoms, 3)

        Returns
        -------
        puw.Quantity
            Scalar
    """
    distance = np.sqrt(np.sum(np.power(coords - centroid, 2), axis=1))
    return np.amax(distance)


def delete(quantity, indices, axis=None):
    """ Return a new array with sub-arrays along an axis deleted.
        For a one dimensional array, this returns those entries not
        returned by arr[indices].

        Parameters
        ----------
        quantity : puw.Quantity

        indices : list[int]

        axis : int, optional

        Returns
        -------
        puw.Quantity
            The new quantity.
    """
    unit = str(puw.get_unit(quantity))
    sub_array = np.delete(puw.get_value(quantity), indices, axis=axis)
    return puw.quantity(sub_array, unit)
