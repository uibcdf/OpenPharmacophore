import numpy as np
from rdkit.Numerics import rdAlignment


def align_pharmacophores(ref_pharma, probe_pharma):
    """ Align two pharmacophores against each other.

        Parameters
        ----------
        ref_pharma : np.ndarray
            Coordinates of the pharmacophoric points of the reference pharmacophore.
            Shape (n_points, 3)

        probe_pharma : np.ndarray
            Coordinates of the pharmacophoric points of probe pharmacophore.
            Shape (n_points, 3)

        Returns
        -------
        rmsd : float
            The root mean square deviation of the alignment.

        trans_mat : np.ndarray
            The transformation matrix. This matrix should be applied to the confomer of
            the probe_pharmacophore to obtain its updated positions.
    """
    ssd, trans_mat = rdAlignment.GetAlignmentTransform(ref_pharma, probe_pharma)
    rmsd = np.sqrt(ssd / ref_pharma.shape[0])
    return rmsd, trans_mat
