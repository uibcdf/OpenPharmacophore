import openpharmacophore as oph
import pyunitwizard as puw
from rdkit import Chem
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D
from rdkit.Chem.Pharm2D.Generate import Gen2DFingerprint


def four_element_pharmacophore():
    pharmacophoric_points = [
        oph.PharmacophoricPoint(
            feat_type="hb acceptor",
            center=puw.quantity([3.877, 7.014, 1.448], "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        ),
        oph.PharmacophoricPoint(
            feat_type="hb acceptor",
            center=puw.quantity([7.22, 11.077, 5.625], "angstroms"),
            radius=puw.quantity(1.0, "angstroms")),
        oph.PharmacophoricPoint(
            feat_type="hb donor",
            center=puw.quantity([4.778, 8.432, 7.805], "angstroms"),
            radius=puw.quantity(1.0, "angstroms")),
        oph.PharmacophoricPoint(
            feat_type="aromatic ring",
            center=puw.quantity([1.56433333333334, 7.06399999999999, 3.135], "angstroms"),
            radius=puw.quantity(1.0, "angstroms"))
    ]
    return oph.Pharmacophore(pharmacophoric_points)


def pharmacophore_fingerprint():
    mol = Chem.MolFromSmiles("Cc1cccc(c2n[nH]cc2c3ccc4ncccc4n3)n1")
    factory = Gobbi_Pharm2D.factory
    return Gen2DFingerprint(mol, factory)
