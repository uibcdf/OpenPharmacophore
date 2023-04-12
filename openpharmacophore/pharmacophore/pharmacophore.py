import numpy as np
import pyunitwizard as puw


class Pharmacophore:
    """ A PharmacophoricPoint container.

        Parameters
        ----------
        points : list[PharmacophoricPoint], optional
            List of pharmacophoric points.

        score : float, optional
            A score assigned to the pharmacophore.

        ref_mol : int, optional
            The index of the reference molecule from which the pharmacophore
            was extracted.

        ref_struct : int, optional
            The index of the reference structure or conformer from which the pharmacophore
            was extracted.

    """
    def __init__(self, points=None, score=None,
                 ref_mol=None, ref_struct=None
                 ):
        if points is None:
            self._points = []
        else:
            self._points = points
        self.score = score
        self.ref_mol = ref_mol
        self.ref_struct = ref_struct
        self.props = {}

    def add(self, point):
        """ Add a pharmacophoric point.

            Parameters
            ----------
            point : PharmacophoricPoint
        """
        self._points.append(point)

    def remove(self, index):
        """ Remove a pharmacophoric point.

            Parameters
            ----------
            index : int
                Index of the pharmacophoric point
        """
        self._points.pop(index)

    def to_matrix(self):
        """ Returns a matrix with the coordinates of the
            pharmacophoric points.
        """
        matrix = np.zeros((len(self), 3))
        for ii in range(len(self)):
            center = puw.get_value(self[ii].center, "angstroms")
            for jj in range(center.shape[0]):
                matrix[ii][jj] = center[jj]

        return puw.quantity(matrix, "angstroms")

    def __eq__(self, other):
        if not self._check_eq(self.score, other.score):
            return False
        if not self._check_eq(self.ref_mol, other.ref_mol):
            return False
        if not self._check_eq(self.ref_struct, other.ref_struct):
            return False

        if len(self) == len(other):
            for ii in range(len(self)):
                if self[ii] != other[ii]:
                    return False
        else:
            return False
        return True

    @staticmethod
    def _check_eq(item1, item2):
        return (item1 is not None and item2 is not None) and item1 == item2

    def __len__(self):
        return len(self._points)

    def __getitem__(self, item):
        return self._points[item]

    def __repr__(self):
        return f"{self.__class__.__name__}(points={self._points}, " \
               f"score={self.score}, ref_mol={self.ref_mol}, ref_struct={self.ref_struct})"
