import openpharmacophore.pharmacophore.chemical_features as feats
import openpharmacophore.data as data
from openpharmacophore import PharmacophoricPoint
from openpharmacophore import utils
import numpy as np
from rdkit import Chem
import pyunitwizard as puw
import pytest


# TODO: Redo this tests

@pytest.fixture
def sample_molecule():
    """Returns a sample molecule for testing"""
    return Chem.MolFromSmiles("CC1=CC(=CC(=C1OCCCC2=CC(=NO2)C)C)C3=NOC(=N3)C(F)(F)F")


def test_feature_centroid(sample_molecule):
    mol = utils.conformers.generate_conformers(molecule=sample_molecule, n_conformers=1, random_seed=1)

    idxs = (11, 12, 13, 14, 15)
    x, y, z = -5.256026015194413, 0.2637169048998521, -0.22204282175815981  # Known coordinates of centroid

    centroid = feats.PharmacophoricPointExtractor._feature_centroid(mol, idxs, 0)

    assert round(x, 2) == round(centroid[0], 2)
    assert round(y, 2) == round(centroid[1], 2)
    assert round(z, 2) == round(centroid[2], 2)


def test_load_smarts_fdef():
    feat_def = feats.load_smarts_fdef(file_name=data.smarts_features)
    assert len(feat_def) == 28

    n_donors = 0
    n_acceptors = 0
    n_positives = 0
    n_negatives = 0
    n_hydrophobes = 0
    for feat_name in feat_def.values():
        if feat_name == "hb donor":
            n_donors += 1
        elif feat_name == "hb acceptor":
            n_acceptors += 1
        elif feat_name == "positive charge":
            n_positives += 1
        elif feat_name == "negative charge":
            n_negatives += 1
        elif feat_name == "hydrophobicity":
            n_hydrophobes += 1

    assert n_donors == 3
    assert n_acceptors == 2
    assert n_hydrophobes == 13
    assert n_negatives == 4
    assert n_positives == 4


def test_pharmacophoric_point_extractor():
    # Known points of acetic acid
    acceptor = PharmacophoricPoint(
        feat_type="hb acceptor",
        center=puw.quantity([0.07903459193908467, 0.1487388014682705, 0.002355596284157232], "nanometer"),
        radius=puw.quantity(1.0, "angstroms")
    )
    negative = PharmacophoricPoint(
        feat_type="negative charge",
        center=puw.quantity([0.08966065176585349, 0.03784658233931392, -0.01482808230401158], "nanometer"),
        radius=puw.quantity(1.0, "angstroms")
    )

    acetic_acid = Chem.MolFromSmiles("CC(=O)O")
    acetic_acid = utils.conformers.generate_conformers(acetic_acid, 1, random_seed=1, alignment=False)

    extractor = feats.PharmacophoricPointExtractor()
    points = extractor(acetic_acid, 0)
    assert len(points) == 3
    assert acceptor.feature_name == points[0].feature_name
    assert np.all(acceptor.center == points[0].center)
    assert acceptor.radius == points[0].radius

    rdkit_extractor = feats.PharmacophoricPointExtractor(featdef=feats.rdkit_featuredefinition())
    points = rdkit_extractor(acetic_acid, 0)
    assert len(points) == 5
    assert acceptor.feature_name == points[0].feature_name
    assert np.all(acceptor.center == points[0].center)
    assert acceptor.radius == points[0].radius

    assert negative.feature_name == points[-1].feature_name
    assert np.all(negative.center == points[-1].center)
    assert negative.radius == points[-1].radius


@pytest.fixture
def benzoic_acid():
    """Returns a benzoic acid molecule for testing"""
    benz_acid = Chem.MolFromSmiles("C1=CC=C(C=C1)C(=O)O")
    benz_acid = utils.conformers.generate_conformers(benz_acid, 1, random_seed=1, alignment=False)
    return benz_acid


def test_aromatic_direction_vector(benzoic_acid):
    aromatic_atoms_inx = (0, 1, 2, 3, 4, 5)
    direction_vector = feats.PharmacophoricPointExtractor._aromatic_direction_vector(
        benzoic_acid, aromatic_atoms_inx, 0)

    assert direction_vector.shape[0] == 3
    assert np.all(np.around(np.array([0.02452289, -1.10048427, 1.20713535]), 2) == np.around(direction_vector, 2))


def test_donor_acceptor_direction_vector(benzoic_acid):
    donor_inx = 8
    donor_center = np.array([2.78976376, 0.87765978, 0.71608437])
    acceptor_inx = 7
    acceptor_center = np.array([2.6851757, -0.79323208, -0.80504238])

    donor_vector = feats.PharmacophoricPointExtractor._donor_acceptor_direction_vector(
        benzoic_acid, 'Donor', donor_inx, donor_center, 0)
    acceptor_vector = feats.PharmacophoricPointExtractor._donor_acceptor_direction_vector(
        benzoic_acid, 'Acceptor', acceptor_inx, acceptor_center, 0)

    assert donor_vector.shape[0] == 3
    assert acceptor_vector.shape[0] == 3
    assert np.all(np.around(np.array([0.7275082, 0.87517266, 0.7830512]), 2) == np.around(donor_vector, 2))
    assert np.all(np.around(np.array([-0.62292014, 0.7957192, 0.73807555]), 2) == np.around(acceptor_vector, 2))
