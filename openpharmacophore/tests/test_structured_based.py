from openpharmacophore import pharmacophoric_elements as phe
from openpharmacophore.structured_based import StructuredBasedPharmacophore as SBP
import pytest
import os

@pytest.mark.parametrize("file_name,ligand_id", [
    ("1ncr", "W11"),
    ("2hz1", None),
    ("2hzi", "JIN")
])
def test_from_pdb(file_name, ligand_id):
    pdbs_path = "./openpharmacophore/data/pdb/"
    file_path = os.path.join(pdbs_path, file_name + ".pdb" )

    pharmacophore = SBP().from_pdb(file_path, radius=1.0, ligand_id=ligand_id)

    if file_name == "1ncr":
       assert len(pharmacophore.elements) == 7
       assert pharmacophore.molecular_system is not None
    elif file_name == "1hz1":
        assert len(pharmacophore.elements) == 12
        assert pharmacophore.molecular_system is not None
    elif file_name == "1hzi":
        assert len(pharmacophore.elements) == 8
        assert pharmacophore.molecular_system is not None
       
@pytest.mark.parametrize("file_name", [
    ("1ncr"),
    ("2hz1"),
    ("2hzi")
])
def test_protein_ligand_interactions(file_name):

    pdbs_path = "./openpharmacophore/data/pdb/"
    file_path = os.path.join(pdbs_path, file_name + ".pdb" )

    interactions, pdb_str = SBP._protein_ligand_interactions(file_path)  
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
    all_interactions, _ = SBP._protein_ligand_interactions(file_path) 

    if file_name == "1ncr":
        interactions = all_interactions["W11:A:7001"]
        pharmacophoric_points = SBP._sb_pharmacophore_points(interactions, radius=1.0)
        assert len(pharmacophoric_points) == 7
        assert isinstance(pharmacophoric_points[0], phe.AromaticRingSphereAndVector)
        n_hydrophobics = 0
        for point in pharmacophoric_points:
            if point.feature_name == "hydrophobicity":
                n_hydrophobics += 1
        assert n_hydrophobics == 6
    
    elif file_name == "2hz1":
        interactions = all_interactions["HEM:A:125"]
        pharmacophoric_points = SBP._sb_pharmacophore_points(interactions, radius=1.0)
        assert len(pharmacophoric_points) == 12
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
        assert n_hydrophobics == 9
        assert n_neg_charges == 2
        assert n_hb_acceptors == 1
    
    elif file_name == "2hzi":
        interactions = all_interactions["JIN:A:600"]
        pharmacophoric_points = SBP._sb_pharmacophore_points(interactions, radius=1.0)
        assert len(pharmacophoric_points) == 8
        assert isinstance(pharmacophoric_points[0], phe.HBAcceptorSphereAndVector)
        assert isinstance(pharmacophoric_points[1], phe.HBDonorSphereAndVector)
        n_hydrophobics = 0
        for point in pharmacophoric_points:
            if point.feature_name == "hydrophobicity":
                n_hydrophobics += 1
        assert n_hydrophobics == 6
    
