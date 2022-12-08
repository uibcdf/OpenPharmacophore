

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

    def __len__(self):
        return len(self._points)

    def __getitem__(self, item):
        return self._points[item]
