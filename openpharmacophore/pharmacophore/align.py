from operator import itemgetter
import numpy as np

from rdkit import Geometry
from rdkit.Chem import rdMolTransforms
from rdkit.Numerics import rdAlignment
from rdkit.Chem.Pharm3D import EmbedLib


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


def align_ligand_to_pharmacophore(ligand, atom_match, pharmacophore):
    """ Align a ligand to a pharmacophore.

        Parameters
        ----------
        ligand : rdkit.Mol
        atom_match : list[list[int]]
        pharmacophore : rdkit.Chem.Pharm3D.Pharmacophore

        Returns
        -------
        rdkit.Mol
            The molecule with an aligned conformer.
        float
            The rmsd of the alignment.
    """

    try:
        _, embeddings, _ = EmbedLib.EmbedPharmacophore(ligand, atom_match, pharmacophore, count=10)
    except KeyError:
        return

    rmsd_list = _transform_embeddings(pharmacophore, embeddings, atom_match)
    if len(rmsd_list) == 0:
        return

    best_fit_index = min(enumerate(rmsd_list), key=itemgetter(1))[0]
    return embeddings[best_fit_index], rmsd_list[best_fit_index]


def _transform_embeddings(pharmacophore, embeddings, atom_match):
    """ Transform embeddings. Performs the alignment of the conformers of
        a molecule to the pharmacophore.

        Parameters
        ----------
        pharmacophore: rdkit.Chem.Pharm3D.Pharmacophore
            A pharmacophore object.

        embeddings: list[rdkit.Mol]
            List of molecules (the same molecule) with a single conformer.

        atom_match: list[list[int]]
            A nested list of atoms ids that match the pharmacophore.

        Returns
        -------
        rmsd_list: list[float]
            List of sum of RMSD values for the alignments.

    """
    align_ref = [f.GetPos() for f in pharmacophore.getFeatures()]
    rmsd_list = []
    for embedding in embeddings:
        conformer = embedding.GetConformer()
        rmsd, transform_matrix = _transform_matrix(align_ref, conformer, atom_match)
        # Transform the coordinates of the conformer
        rdMolTransforms.TransformConformer(conformer, transform_matrix)
        rmsd_list.append(rmsd)

    return rmsd_list


def _transform_matrix(align_ref, conformer_embed, atom_match):
    """ Get the transformation matrix of a conformer that is aligned to
        a pharmacophore.

        Parameters
        ----------
        align_ref: list[rdkit.Geometry.Point3D]
            list of pharmacophore reference points for the alignment.

        conformer_embed: rdkit.Conformer
            The conformer embedding.

        atom_match: list[list]
            Nested list of atoms ids that match the pharmacophore.

        Returns
        -------
        rmsd: float
            RMSD value of the alignment.

        transform_matrix: numpy.ndarray; shape(4, 4)
            The transform matrix.
    """
    align_probe = []  # probe points to align to reference points
    for match_ids in atom_match:
        # Calculate the centroid of the feature in case it has multiple atoms
        dummy_point = Geometry.Point3D(0.0, 0.0, 0.0)
        for id_ in match_ids:
            dummy_point += conformer_embed.GetAtomPosition(id_)
        dummy_point /= len(match_ids)
        align_probe.append(dummy_point)

    # Note This function returns the sum of squared distance (SSD) not the RMSD RMSD = sqrt(SSD/numPoints)
    ssd, transform_matrix = rdAlignment.GetAlignmentTransform(align_ref, align_probe)
    rmsd = np.sqrt(ssd / len(align_ref))
    return rmsd, transform_matrix
