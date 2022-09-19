import pyunitwizard as puw
from rdkit import Geometry
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm3D import Pharmacophore as rdkitPharmacophore


def rdkit_pharmacophore(pharmacophoric_points):
    """ Transform a list of pharmacophoric points to a rdkit pharmacophore.

        rdkit pharmacophores do not store the pharmacophoric_points radii, so they are returned as well.

        Returns
        -------
        rdkit_pharmacophore : rdkit.Chem.Pharm3D.Pharmacophore
            The rdkit pharmacophore.

        radii : list[float]
            List with the radius in angstroms of each pharmacophoric point.
    """
    rdkit_element_name = {
        "aromatic ring": "Aromatic",
        "hydrophobicity": "Hydrophobe",
        "hb acceptor": "Acceptor",
        "hb donor": "Donor",
        "positive charge": "PosIonizable",
        "negative charge": "NegIonizable",
    }

    points = []
    radii = []

    for element in pharmacophoric_points:
        feat_name = rdkit_element_name[element.feature_name]
        center = puw.get_value(element.center, to_unit="angstroms")
        center = Geometry.Point3D(center[0], center[1], center[2])
        points.append(ChemicalFeatures.FreeChemicalFeature(
            feat_name,
            center
        ))
        radius = puw.get_value(element.radius, to_unit="angstroms")
        radii.append(radius)

    pharmacophore = rdkitPharmacophore.Pharmacophore(points)
    # Apply the radius of each point to the bound matrix of the pharmacophore
    for ii in range(len(radii)):
        for jj in range(ii + 1, len(radii)):
            sum_radii = radii[ii] + radii[jj]
            pharmacophore.setLowerBound(ii, jj, max(pharmacophore.getLowerBound(ii, jj) - sum_radii, 0))
            pharmacophore.setUpperBound(ii, jj, pharmacophore.getUpperBound(ii, jj) + sum_radii)

    return pharmacophore
