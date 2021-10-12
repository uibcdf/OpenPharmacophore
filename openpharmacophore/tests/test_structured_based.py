from openpharmacophore.pharmacophoric_point import PharmacophoricPoint
from openpharmacophore.structured_based import StructuredBasedPharmacophore as SBP
from openpharmacophore.utils.conformers import generate_conformers
import numpy as np
import pytest
import pyunitwizard as puw
from rdkit import Chem
import os

from openpharmacophore.utils.conformers import generate_conformers

@pytest.mark.parametrize("file_name,ligand_id,hydrophobics", [
    ("1ncr", "W11:A:7001", "plip"),
    ("2hz1", None, "plip"),
    ("2hzi", "JIN:A:600", "plip"),
    ("1qku", "EST:A:600", "rdkit")
])
def test_from_pdb(file_name, ligand_id, hydrophobics):
    pdbs_path = "./openpharmacophore/data/pdb/"
    file_path = os.path.join(pdbs_path, file_name + ".pdb" )

    pharmacophore = SBP().from_pdb(file_path, radius=1.0, ligand_id=ligand_id, hydrophobics=hydrophobics)

    if file_name == "1ncr":
        assert len(pharmacophore.elements) == 5
        assert pharmacophore.molecular_system is not None
        
        ring = PharmacophoricPoint(
            feat_type="aromatic ring",
            center=puw.quantity((38.0292, 4.3594, 123.7382), "angstroms"),
            radius=puw.quantity(1.0, "angstroms"),
            direction = np.array([-0.6702, -0.0512, -0.7404])
        )
        assert ring == pharmacophore.elements[0]
        hyd_1 = PharmacophoricPoint(
            feat_type="hydrophobicity",
            center=puw.quantity((42.4462, -1.1092, 122.6226), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert hyd_1 == pharmacophore.elements[1]
        hyd_2 = PharmacophoricPoint(
            feat_type="hydrophobicity",
            center=puw.quantity((37.699, 4.393, 122.609), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert hyd_2 == pharmacophore.elements[2]
        hyd_3 = PharmacophoricPoint(
            feat_type="hydrophobicity",
            center=puw.quantity((35.827, 5.861, 123.815), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert hyd_3 == pharmacophore.elements[3]
        hyd_4 = PharmacophoricPoint(
            feat_type="hydrophobicity",
            center=puw.quantity((45.46, 2.006, 123.421), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert hyd_4 == pharmacophore.elements[4]
    elif file_name == "2hz1":
        assert len(pharmacophore.elements) == 10
        assert pharmacophore.molecular_system is not None
        assert pharmacophore.ligand is None
        neg_sphere = PharmacophoricPoint(
            feat_type="negative charge",
            center=puw.quantity((2.279, 2.6625, 3.2525), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert neg_sphere == pharmacophore.elements[0]
        neg_sphere = PharmacophoricPoint(
            feat_type="negative charge",
            center=puw.quantity((1.074, 4.8205, 7.3245), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert neg_sphere == pharmacophore.elements[1]
        hb_acceptor = PharmacophoricPoint(
            feat_type="hb acceptor",
            center=puw.quantity((2.327, 3.574, 2.628), "angstroms"),
            radius=puw.quantity(1.0, "angstroms"),
            direction = np.array([-0.8834, 0.2719, 0.3816])
        )
        assert hb_acceptor == pharmacophore.elements[2]
        hyd = PharmacophoricPoint(
            feat_type="hydrophobicity",
            center=puw.quantity((5.134, 11.9745, 11.8155), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert hyd == pharmacophore.elements[3]
        hyd = PharmacophoricPoint(
            feat_type="hydrophobicity",
            center=puw.quantity((0.396, 6.793, 7.373), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert hyd == pharmacophore.elements[4]
        hyd = PharmacophoricPoint(
            feat_type="hydrophobicity",
            center=puw.quantity((11.724, 5.702, 7.226), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert hyd == pharmacophore.elements[5]
        hyd = PharmacophoricPoint(
            feat_type="hydrophobicity",
            center=puw.quantity((12.502, 7.731, 9.51), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert hyd == pharmacophore.elements[6]
        hyd = PharmacophoricPoint(
            feat_type="hydrophobicity",
            center=puw.quantity((7.62, 3.602, 4.389), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert hyd == pharmacophore.elements[7]
        hyd = PharmacophoricPoint(
            feat_type="hydrophobicity",
            center=puw.quantity((8.728, 11.736, 11.502), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert hyd == pharmacophore.elements[8]

    elif file_name == "2hzi":
        assert len(pharmacophore.elements) == 5
        assert pharmacophore.molecular_system is not None
        assert pharmacophore.ligand is None
        acceptor = PharmacophoricPoint(
            feat_type="hb acceptor",
            center=puw.quantity((16.163, 12.126, 2.521), "angstroms"),
            radius=puw.quantity(1.0, "angstroms"),
            direction = np.array([-0.0991, 0.9107, 0.4011])
        )
        assert acceptor == pharmacophore.elements[0]
        donor = PharmacophoricPoint(
            feat_type="hb donor",
            center=puw.quantity((13.988, 12.055, 2.867), "angstroms"),
            radius=puw.quantity(1.0, "angstroms"),
            direction = np.array([-0.1195, -0.9419, -0.3139])
        )
        assert donor == pharmacophore.elements[1]
        hyd = PharmacophoricPoint(
            feat_type="hydrophobicity",
            center=puw.quantity((10.418, 11.8007, 2.3763), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert hyd == pharmacophore.elements[2]
        hyd = PharmacophoricPoint(
            feat_type="hydrophobicity",
            center=puw.quantity((21.7962, 17.734, 4.249), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert hyd == pharmacophore.elements[3]
        hyd = PharmacophoricPoint(
            feat_type="hydrophobicity",
            center=puw.quantity((17.506, 13.853, 3.414), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert hyd == pharmacophore.elements[4]
    elif file_name == "1qku":
        assert len(pharmacophore.elements) == 6
        assert pharmacophore.molecular_system is not None
        assert pharmacophore.ligand is not None
        acceptor = PharmacophoricPoint(
            feat_type="hb acceptor",
            center=puw.quantity((100.598, 17.99, 25.704), "angstroms"),
            radius=puw.quantity(1.0, "angstroms"),
            direction = np.array([0.7872, -0.4533, 0.4182])
        )
        assert acceptor == pharmacophore.elements[0]
        acceptor = PharmacophoricPoint(
            feat_type="hb acceptor",
            center=puw.quantity((109.61, 12.027, 22.746), "angstroms"),
            radius=puw.quantity(1.0, "angstroms"),
            direction = np.array([-0.5805, 0.431, 0.6908])
        )
        assert acceptor == pharmacophore.elements[1]
        ring = PharmacophoricPoint(
            feat_type="aromatic ring",
            center=puw.quantity((102.8133, 16.7163, 24.5195), "angstroms"),
            radius=puw.quantity(1.0, "angstroms"),
            direction = np.array([0.0937, 0.4983, -0.8619])
        )
        assert ring == pharmacophore.elements[3]
        hyd = PharmacophoricPoint(
            feat_type="hydrophobicity",
            center=puw.quantity((107.379, 12.449, 24.419), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert hyd == pharmacophore.elements[4]
        hyd = PharmacophoricPoint(
            feat_type="hydrophobicity",
            center=puw.quantity((107.2112, 12.5732, 22.1112), "angstroms"),
            radius=puw.quantity(1.0, "angstroms")
        )
        assert hyd == pharmacophore.elements[5]

       
@pytest.mark.parametrize("file_name", [
    ("1ncr"),
    ("2hz1"),
    ("2hzi")
])
def test_protein_ligand_interactions(file_name):

    pdbs_path = "./openpharmacophore/data/pdb/"
    file_path = os.path.join(pdbs_path, file_name + ".pdb" )

    interactions, pdb_str, _ = SBP._protein_ligand_interactions(file_path, as_string=False)  
    assert isinstance(interactions, dict)
    assert isinstance(pdb_str, str)
    
    ligands = list(interactions.keys())

    if file_name == "1ncr":
        assert len(interactions) == 2
        assert ligands[0] == "W11:A:7001"
        assert ligands[1] == "MYR:D:4000"
    elif file_name == "1hz1":
        assert len(interactions) == 1
        assert ligands[0] == "HEM:A:125"
    elif file_name == "1hzi":
        assert len(interactions) == 2
        assert ligands[0] == "JIN:A:600"
        assert ligands[1] == "JIN:B:600"

@pytest.mark.parametrize("file_name", [
    ("1ncr"),
    ("2hz1"),
    ("2hzi")
])
def test_sb_pharmacophore_points(file_name):
    pdbs_path = "./openpharmacophore/data/pdb/"
    file_path = os.path.join(pdbs_path, file_name + ".pdb" )
    all_interactions, _ , _ = SBP._protein_ligand_interactions(file_path, as_string=False) 

    if file_name == "1ncr":
        interactions = all_interactions["W11:A:7001"]
        pharmacophoric_points = SBP._sb_pharmacophore_points(interactions, radius=1.0, ligand=None, hydrophobics="plip")
        assert len(pharmacophoric_points) == 5
        assert isinstance(pharmacophoric_points[0], PharmacophoricPoint)
        assert pharmacophoric_points[0].feature_name == "aromatic ring"
        n_hydrophobics = 0
        for point in pharmacophoric_points:
            if point.feature_name == "hydrophobicity":
                n_hydrophobics += 1
        assert n_hydrophobics == 4
    
    elif file_name == "2hz1":
        interactions = all_interactions["HEM:A:125"]
        pharmacophoric_points = SBP._sb_pharmacophore_points(interactions, radius=1.0, ligand=None, hydrophobics="plip")
        assert len(pharmacophoric_points) == 10
        n_hydrophobics = 0
        n_neg_charges = 0
        n_hb_acceptors = 0
        for point in pharmacophoric_points:
            if point.feature_name == "hydrophobicity":
                n_hydrophobics += 1
            elif point.feature_name == "hb acceptor":
                n_hb_acceptors += 1
            elif point.feature_name == "negative charge":
                n_neg_charges += 1
        assert n_hydrophobics == 7
        assert n_neg_charges == 2
        assert n_hb_acceptors == 1
    
    elif file_name == "2hzi":
        interactions = all_interactions["JIN:A:600"]
        pharmacophoric_points = SBP._sb_pharmacophore_points(interactions, radius=1.0, ligand=None, hydrophobics="plip")
        assert len(pharmacophoric_points) == 5
        assert isinstance(pharmacophoric_points[0], PharmacophoricPoint)
        assert pharmacophoric_points[0].feature_name == "hb acceptor"
        assert isinstance(pharmacophoric_points[1], PharmacophoricPoint)
        assert pharmacophoric_points[1].feature_name == "hb donor"
        n_hydrophobics = 0
        for point in pharmacophoric_points:
            if point.feature_name == "hydrophobicity":
                n_hydrophobics += 1
        assert n_hydrophobics == 3
    
def test_rdkit_hydrophobics():
    
    ligand = Chem.MolFromSmiles("CC1=CC(=CC(=C1OCCCC2=CC(=NO2)C)C)C3=NOC(=N3)C(F)(F)F")
    ligand = generate_conformers(ligand, 1, random_seed=1)
    hydrophobics = SBP()._rdkit_hydrophobics(ligand, 1.0)

    assert len(hydrophobics) == 2
    
    hyd_1 = hydrophobics[0]
    hyd_2 = hydrophobics[1]
    hyd_2_center = puw.get_value(hyd_2.center, "angstroms")
    hyd_2_radius = puw.get_value(hyd_2.radius, "angstroms")
    assert hyd_1 == PharmacophoricPoint(
                    feat_type="hydrophobicity",
                    center=puw.quantity((1.4268, -1.9841, -1.4186), "angstroms"),
                    radius=puw.quantity(1.0, "angstroms")
                    )
    assert np.allclose(np.around(hyd_2_center, 3), (-2.130, -0.542, -0.479), rtol=0, atol=1e-04)
    assert np.allclose(hyd_2_radius, hyd_2_radius, rtol=0, atol=1e-03) 