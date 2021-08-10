from openpharmacophore.utils.load_custom_feats import load_smarts_fdef
from openpharmacophore._private_tools.exceptions import PointTypeError
from openpharmacophore.utils.rdkit_to_point import rdkit_to_point
from openpharmacophore import pharmacophoric_elements, utils
from openpharmacophore.utils.ligand_features import ligands_pharmacophoric_points
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

    # Round to two decimals because precission varies from machine to machine
    assert round(x, 2) == round(centroid[0], 2)
    assert round(y, 2) == round(centroid[1], 2)
    assert round(z, 2) == round(centroid[2], 2)


@pytest.mark.parametrize("radius,feat_def,point_type,exception",[
    (1.0, "rdkit", "spheres", None),
    (1.0, "rdkit", "ball", PointTypeError("Invalid point type. \"ball\" is not a valid point type")),
    (None, "rdkit", "spheres_vectors", ValueError("Radius cannot be null if point type is spheres or spheres_vectors")),
    (1.0, "custom", "spheres", None)
])
def test_ligands_pharmacophoric_points(radius, feat_def, point_type, exception):
    # Known points of acetic acid
    acceptor = pharmacophoric_elements.HBAcceptorSphere(
        center=puw.quantity([0.07903459193908467, 0.1487388014682705, 0.002355596284157232], "nanometer"),
        radius=puw.quantity(1.0, "angstroms")
    )
    negative = pharmacophoric_elements.NegativeChargeSphere(
        center=puw.quantity([0.08966065176585349, 0.03784658233931392, -0.01482808230401158], "nanometer"),
        radius=puw.quantity(1.0, "angstroms")
    )
    
    acetic_acid = Chem.MolFromSmiles("CC(=O)O")
    acetic_acid = utils.conformers.generate_conformers(acetic_acid, 1, random_seed=1, alignment=False)
    
    if feat_def == "rdkit":
        definition = None
    elif feat_def == "custom":
        definition = load_smarts_fdef(fname="openpharmacophore/data/smarts_features.txt")
    
    try:
        points = ligands_pharmacophoric_points(acetic_acid, feat_list=None, feat_def=definition, 
                            point_type=point_type, radius=radius)
    except (PointTypeError, ValueError) as e:
        e.args == exception.args
    else:
        points = points["ligand_0"]["conformer_0"]

        if feat_def == "rdkit":
            assert len(points) == 5 
            assert acceptor.feature_name == points[1].feature_name
            assert np.all(acceptor.center == points[1].center)
            assert acceptor.radius == points[1].radius
            assert acceptor.shape_name == points[1].shape_name

            assert negative.feature_name == points[3].feature_name
            assert np.all(negative.center == points[3].center)
            assert negative.radius == points[3].radius
            assert negative.shape_name == points[3].shape_name
        
        elif feat_def == "custom":
            assert len(points) == 3
            assert acceptor.feature_name == points[0].feature_name
            assert np.all(acceptor.center == points[0].center)
            assert acceptor.radius == points[0].radius
            assert acceptor.shape_name == points[0].shape_name

def test_rdkit_to_point():
    aromatic_sphere = rdkit_to_point("Aromatic", [0.0, 1.0, 1.0], radius=1.0, point_type="spheres")  
    assert isinstance(aromatic_sphere, pharmacophoric_elements.AromaticRingSphere)
    assert aromatic_sphere.radius == puw.standardize(puw.quantity(1.0, "angstroms"))
    assert aromatic_sphere.feature_name == "aromatic ring"
    assert np.all(aromatic_sphere.center == puw.standardize(puw.quantity([0.0, 1.0, 1.0], "angstroms")))

    donor_vector = rdkit_to_point("Donor", [-1.0, 0.0, -1.0], radius=1.0, direction=[1.0, 1.0, 1.0], point_type="spheres_vectors")
    assert isinstance(donor_vector, pharmacophoric_elements.HBDonorSphereAndVector)
    assert donor_vector.radius == puw.standardize(puw.quantity(1, "angstroms"))
    assert donor_vector.feature_name == "hb donor"
    assert np.all(donor_vector.center == puw.standardize(puw.quantity([-1.0, 0.0, -1.0], "angstroms")))
    assert np.all(donor_vector.direction == (np.array([1.0, 1.0, 1.0]) / np.linalg.norm([1.0, 1.0, 1.0])))

    hydrophobe_gaussian = rdkit_to_point("Hydrophobe", [1.0, 0.0, 0.0], sigma=0.8, point_type="gaussian")
    assert isinstance(hydrophobe_gaussian, pharmacophoric_elements.HydrophobicGaussianKernel)
    assert hydrophobe_gaussian.sigma == puw.standardize(puw.quantity(0.8, "angstroms"))
    assert hydrophobe_gaussian.feature_name == "hydrophobicity"
    assert np.all(hydrophobe_gaussian.center == puw.standardize(puw.quantity([1.0, 0.0, 0.0], "angstroms")))


def test_load_smarts_fdef():
    feat_def = load_smarts_fdef(fname="openpharmacophore/data/smarts_features.txt")

    assert len(feat_def) == 28
    assert '[#16!H0]' in feat_def
    assert 'c1nn[nH1]n1' in feat_def

def test_direction_vector():
    dir_vector = [1, 0, 0]
    assert np.all(dir_vector == direction_vector())