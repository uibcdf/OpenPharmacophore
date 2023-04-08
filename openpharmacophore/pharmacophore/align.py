import numpy as np


def align_pharmacophores(pharma1, pharma2):
    """ Align two pharmacophores against each other.

        Parameters
        ----------
        pharma1 : QuantityLike
            Coordinates of the pharmacophoric points of first pharmacophore.

        pharma2 : QuantityLike
            Coordinates of the pharmacophoric points of second pharmacophore.

        Returns
        -------
        rmsd : float
            The root mean square deviation of the alignment.
    """
    # TODO: implement algorithm to align pharmacophores
    rmsd = np.sqrt(np.power(pharma1 - pharma2, 2).mean())
    return rmsd
