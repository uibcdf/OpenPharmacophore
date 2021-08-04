from openpharmacophore.utils.rdkit_to_point import rdkit_to_point
from openpharmacophore import pharmacophoric_elements, utils
from openpharmacophore.utils.ligand_features import ligands_pharmacophoric_points
from rdkit import Chem
import pyunitwizard as puw
import numpy as np

def test_generate_conformers():
    mol = Chem.MolFromSmiles("CC1=CC(=CC(=C1OCCCC2=CC(=NO2)C)C)C3=NOC(=N3)C(F)(F)F")
    mol = utils.conformers.generate_conformers(molecule=mol, n_conformers=4)
    assert mol.GetNumConformers() == 4

# def test_feature_centroid():
#     mol = Chem.MolFromSmiles("CC1=CC(=CC(=C1OCCCC2=CC(=NO2)C)C)C3=NOC(=N3)C(F)(F)F")
#     mol = utils.conformers.generate_conformers(molecule=mol, n_conformers=1, random_seed=1)
# 
#     idxs = (11, 12, 13, 14, 15)
#     x, y, z = -5.256026015194413, 0.2637169048998521, -0.22204282175815981 # Known coordinates of centroid
# 
#     assert x == utils.centroid.feature_centroid(mol, idxs, 0)[0]
#     assert y == utils.centroid.feature_centroid(mol, idxs, 0)[1]
#     assert z == utils.centroid.feature_centroid(mol, idxs, 0)[2]

def test_ligands_pharmacophoric_points():
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
    points = ligands_pharmacophoric_points(acetic_acid, feat_list=None, feat_def='rdkit', 
                               point_type="spheres", radius=1)
    

    points = points["ligand_0"]["conformer_0"]

    assert acceptor.feature_name == points[0].feature_name
    assert np.all(acceptor.center == points[0].center)
    assert acceptor.radius == points[0].radius
    assert acceptor.shape_name == points[0].shape_name

    assert negative.feature_name == points[2].feature_name
    assert np.all(negative.center == points[2].center)
    assert negative.radius == points[2].radius
    assert negative.shape_name == points[2].shape_name

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
