# This file contains functions necessary to align molecules to a 3D pharmacophore
# model. The functions in this file are called by the VirtualScreening3D and
# RetrospectiveScreening3D classes

from rdkit import Geometry
from rdkit.Chem import rdMolTransforms
from rdkit.Numerics import rdAlignment


def apply_radii_to_bounds(radii, pharmacophore):
    """
        Apply the radius of the pharmacophoric points to the bound
        matrix of the pharmacophore

        Parameters
        ----------
        radii: list of float
            list with the radius of each of the pharmacophoric points
        
        pharmacophore: rdkit.Chem.Pharm3D.Pharmacophore
            The pharmacophore wich bounds matriz will be modified
        
        Returns
        -------
        Nothing is returned. The pharmacophore bounds matrix is updated

    """
    for i in range(len(radii)):
        for j in range(i + 1, len(radii)):
            sum_radii = radii[i] + radii[j]
            pharmacophore.setLowerBound(i, j, max(pharmacophore.getLowerBound(i, j) - sum_radii, 0))
            pharmacophore.setUpperBound(i, j, pharmacophore.getUpperBound(i, j) + sum_radii)


def get_transform_matrix(align_ref, conf_embed, atom_match):

    """Get the transformation matrix for a conformer.

        Parameters
        ----------
        align_ref: list of rdkit.Geometry.Point3D
            list of pharmacophore reference points for the alignment.

        conf_embed: rdkit.Chem.Conformer
            The conformer embedding.

        atom_match: list of list
            List of list of atoms ids that match the pharmacophore.

        Returns
        -------
        a 2-tuple

        ssd: float
            SSD value for the alignment.

        transform_matrix: numpy.ndarray; shape(4, 4)
            The transform matrix.

        """

    align_probe = [] # probe points to align to reference points
    for match_ids in atom_match:
        dummy_point = Geometry.Point3D(0.0,0.0,0.0)
        for id in match_ids:
            dummy_point += conf_embed.GetAtomPosition(id)
        dummy_point /= len(match_ids)
        align_probe.append(dummy_point)
    
    ssd, transform_matrix = rdAlignment.GetAlignmentTransform(align_ref, align_probe)
    return ssd, transform_matrix


def transform_embeddings(pharmacophore, embeddings, atom_match):
    """Transform embeddings. Performs the alignment of the molecules 
        to the pharmacophore.

        Parameters
        ----------
        pharmacophore: rdkit.Chem.Pharm3D.Pharmacophore
            A pharmacophore object.

        embeddings: list[rdkit.Mol]
            List of molecules with a single conformer.

        atom_match: list[list[int]]
            A nested list of atoms ids that match the pharmacophore.

        Returns
        -------
        SSDs: list[float]
            List of sum of square deviations (SSD) values for the alignments.

        """

    align_ref = [f.GetPos() for f in pharmacophore.getFeatures()]
    ssds = []
    for embedding in embeddings:
        conformer = embedding.GetConformer()
        ssd, transform_matrix = get_transform_matrix(align_ref, conformer, atom_match)
        # Transform the coordinates of the conformer
        rdMolTransforms.TransformConformer(conformer, transform_matrix)
        ssds.append(ssd)
    
    return ssds
