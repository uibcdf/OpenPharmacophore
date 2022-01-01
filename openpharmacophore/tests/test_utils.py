from openpharmacophore.pharmacophoric_point import PharmacophoricPoint
from openpharmacophore.utils.direction_vector import aromatic_direction_vector, donor_acceptor_direction_vector
from openpharmacophore.utils.load_custom_feats import load_smarts_fdef
from openpharmacophore.algorithms.discretize import discretize
from openpharmacophore import utils
from openpharmacophore.utils.ligand_features import ligands_pharmacophoric_points, rdkit_to_point
from rdkit import Chem
import pyunitwizard as puw
import numpy as np
import pytest

@pytest.fixture
def sample_molecule():
    """Returns a sample molecule for testing"""
    return Chem.MolFromSmiles("CC1=CC(=CC(=C1OCCCC2=CC(=NO2)C)C)C3=NOC(=N3)C(F)(F)F")

def test_generate_conformers(sample_molecule):
    mol = utils.conformers.generate_conformers(molecule=sample_molecule, n_conformers=2)
    assert mol.GetNumConformers() == 2

def test_feature_centroid(sample_molecule):
    mol = utils.conformers.generate_conformers(molecule=sample_molecule, n_conformers=1, random_seed=1)

    idxs = (11, 12, 13, 14, 15)
    x, y, z = -5.256026015194413, 0.2637169048998521, -0.22204282175815981 # Known coordinates of centroid

    centroid = utils.centroid.feature_centroid(mol, idxs, 0)

    assert round(x, 2) == round(centroid[0], 2)
    assert round(y, 2) == round(centroid[1], 2)
    assert round(z, 2) == round(centroid[2], 2)


@pytest.mark.parametrize("radius,feat_def",[
    (1.0, "rdkit"),
    (1.0, "custom")
])
def test_ligands_pharmacophoric_points(radius, feat_def,):
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
    
    if feat_def == "rdkit":
        definition = None
    elif feat_def == "custom":
        definition = load_smarts_fdef(file_name="openpharmacophore/data/smarts_features.txt")
    
    points = ligands_pharmacophoric_points(
                    acetic_acid, 
                    feat_list=None,
                    feat_def=definition, 
                    radius=radius)

    points = points["ligand_0"]["conformer_0"]

    if feat_def == "rdkit":
        assert len(points) == 5 
        assert acceptor.feature_name == points[1].feature_name
        assert np.all(acceptor.center == points[1].center)
        assert acceptor.radius == points[1].radius

        assert negative.feature_name == points[3].feature_name
        assert np.all(negative.center == points[3].center)
        assert negative.radius == points[3].radius
    
    elif feat_def == "custom":
        assert len(points) == 3
        assert acceptor.feature_name == points[0].feature_name
        assert np.all(acceptor.center == points[0].center)
        assert acceptor.radius == points[0].radius

def test_rdkit_to_point():
    aromatic_sphere = rdkit_to_point("Aromatic", [0.0, 1.0, 1.0], radius=1.0)  
    assert isinstance(aromatic_sphere, PharmacophoricPoint)
    assert aromatic_sphere.radius == puw.standardize(puw.quantity(1.0, "angstroms"))
    assert aromatic_sphere.feature_name == "aromatic ring"
    assert np.all(aromatic_sphere.center == puw.standardize(puw.quantity([0.0, 1.0, 1.0], "angstroms")))

    donor_vector = rdkit_to_point("Donor", [-1.0, 0.0, -1.0], radius=1.0, direction=[1.0, 1.0, 1.0])
    assert isinstance(donor_vector, PharmacophoricPoint)
    assert donor_vector.radius == puw.standardize(puw.quantity(1, "angstroms"))
    assert donor_vector.feature_name == "hb donor"
    assert np.all(donor_vector.center == puw.standardize(puw.quantity([-1.0, 0.0, -1.0], "angstroms")))
    assert np.all(donor_vector.direction == (np.array([1.0, 1.0, 1.0]) / np.linalg.norm([1.0, 1.0, 1.0])))


def test_load_smarts_fdef():
    feat_def = load_smarts_fdef(file_name="openpharmacophore/data/smarts_features.txt")

    assert len(feat_def) == 28
    assert '[#16!H0]' in feat_def
    assert 'c1nn[nH1]n1' in feat_def

@pytest.fixture
def benzoic_acid():
    """Returns a benzoic acid molecule for testing"""
    benz_acid = Chem.MolFromSmiles("C1=CC=C(C=C1)C(=O)O")
    benz_acid = utils.conformers.generate_conformers(benz_acid, 1, random_seed=1, alignment=False)
    return benz_acid

def test_aromatic_direction_vector(benzoic_acid):

    aromatic_atoms_inx = (0, 1, 2, 3, 4, 5)
    direction_vector = aromatic_direction_vector(benzoic_acid, aromatic_atoms_inx, 0)

    assert direction_vector.shape[0] == 3
    assert np.all(np.around(np.array([ 0.02452289, -1.10048427,  1.20713535]), 2) == np.around(direction_vector, 2))

def test_donor_acceptor_direction_vector(benzoic_acid):
    donor_inx = 8
    donor_center = np.array([2.78976376, 0.87765978, 0.71608437])
    acceptor_inx = 7
    acceptor_center = np.array([2.6851757, -0.79323208, -0.80504238])

    donor_vector = donor_acceptor_direction_vector(benzoic_acid, 'Donor', donor_inx, donor_center, 0)
    acceptor_vector = donor_acceptor_direction_vector(benzoic_acid, 'Acceptor', acceptor_inx, acceptor_center, 0)

    assert donor_vector.shape[0] == 3
    assert acceptor_vector.shape[0] == 3
    assert np.all(np.around(np.array([0.7275082, 0.87517266, 0.7830512]), 2) == np.around(donor_vector, 2))
    assert np.all(np.around(np.array([-0.62292014, 0.7957192 , 0.73807555]), 2) == np.around(acceptor_vector, 2))
    
    
def test_discretize():
    bins = np.arange(0, 10, 1.0)
    
    assert discretize(5.99999, bins) == 6.0
    assert discretize(6.0, bins) == 6.0
    assert discretize(5.5, bins) == 5.0
    assert discretize(5.4, bins) == 5.0
    assert discretize(5.6, bins) == 6.0
    