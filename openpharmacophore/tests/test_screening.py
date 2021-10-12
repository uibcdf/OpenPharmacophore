from openpharmacophore.pharmacophore import Pharmacophore
from openpharmacophore.pharmacophoric_point import PharmacophoricPoint
from openpharmacophore.screening import screening, screening2D, screening3D
import numpy as np
import pytest
import pyunitwizard as puw
import pandas as pd
from rdkit import Chem
import os

### Tests for VirtrualScreening base class ###

@pytest.fixture
def mock_screening_results():
    """ Returns an isntance of a VirtualScreening class
        with fake attributes.
    """
    pharmacophore = Pharmacophore()
    screener = screening.VirtualScreening(pharmacophore)
    screener.n_molecules = 60000
    screener.n_matches = 3
    screener.n_fails = 60000 - 3
    screener.scoring_metric = "Similarity"
    screener.db = "ChemBl"

    screener.matches = [
        (0.6542, "CHEMBL22796", Chem.MolFromSmiles("CCCCc1nn(CCC)c(C(=O)O)c1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1")),
        (0.7833, "CHEMBL431887", Chem.MolFromSmiles("CCCCc1nn(CC)c(C(=O)O)c1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1")),
        (0.8974, "CHEMBL907", Chem.MolFromSmiles("CCCCc1nc(Cl)c(C(=O)O)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1")),
    ]

    return screener

def test_get_report(mock_screening_results):
    screener = mock_screening_results
    report_str = screener._get_report()
    print(report_str)

    expected_report = """
        Virtual Screening Results
        -------------------------

        Molecules scanned:                               60,000
        Molecules matched to pharmacophore:                   3
        Molecules that didn't match the pharmacophore:   59,997
        Lowest  Similarity value:     0.6542
        Highest Similarity value:     0.8974
        Average Similarity value:     0.7783

        Top 3 molecules:

        ChemBlID     Similarity
        -------       ------
        CHEMBL907     0.8974
        CHEMBL431887  0.7833
        CHEMBL22796   0.6542""".split()

    assert expected_report == report_str.split()

@pytest.mark.parametrize("form", ["dict", "dataframe"])
def test_get__screening_results(form, mock_screening_results):
    screener = mock_screening_results
    results = screener.get_screening_results(form=form)

    if form == "dict":
        assert isinstance(results, dict)
        assert len(results) == 5
        assert "ChemBl_id" in results
        assert results["ChemBl_id"] == ["CHEMBL22796",
                            "CHEMBL431887", "CHEMBL907"]
        assert results["Smiles"] == ["CCCCc1nn(CCC)c(C(=O)O)c1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1",
            "CCCCc1nn(CC)c(C(=O)O)c1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1",
            "CCCCc1nc(Cl)c(C(=O)O)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1"
        ]
        assert np.allclose(np.array(results["Mol_weight"]), np.array([444.54, 430.51, 436.90]))
        assert len(results["logP"]) == 3
        assert results["Similarity"] == [0.6542, 0.7833, 0.8974]

    else:
        assert isinstance(results, pd.DataFrame)
        assert list(results.columns) == ['ChemBl_id', 'Smiles', 'Similarity', 'Mol_weight', 'logP']
        assert results["ChemBl_id"].to_list() == ["CHEMBL22796",  "CHEMBL431887", "CHEMBL907"]
        assert results["Smiles"].to_list() == ["CCCCc1nn(CCC)c(C(=O)O)c1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1",
            "CCCCc1nn(CC)c(C(=O)O)c1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1",
            "CCCCc1nc(Cl)c(C(=O)O)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1"
        ]
        assert np.allclose(results["Mol_weight"].to_numpy(), np.array([444.54, 430.51, 436.90]))
        assert np.allclose(results["logP"].to_numpy(), np.array([4.7718, 4.3817, 4.4727]))
        assert np.allclose(results["Similarity"].to_numpy(), np.array([0.6542, 0.7833, 0.8974]))

@pytest.mark.parametrize("file_name", [
    "ace.mol2",
    "clique_detection.smi"
])
def test_load_molecules_file(file_name):
    pharmacophore = Pharmacophore()
    screener = screening.VirtualScreening(pharmacophore)
    file_path = os.path.join("./openpharmacophore/data/ligands", file_name)
    
    if file_name.endswith(".smi"):
        ligands = screener._load_molecules_file(file_name=file_path, titleLine=False)
        assert len(ligands) == 5
    elif file_name.endswith(".mol2"):
        ligands = screener._load_molecules_file(file_name=file_path)
        assert len(ligands) == 3

    assert isinstance(ligands, list)
    for lig in ligands:
        assert isinstance(lig, Chem.Mol)

### Tests for VirtrualScreening3D class ###
def test_screen_db_from_dir_3D():

    elements = [
        PharmacophoricPoint(
        feat_type="hb acceptor",
        center=puw.quantity([3.877, 7.014, 1.448], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
        ),
        PharmacophoricPoint(
        feat_type="hb acceptor",
        center=puw.quantity([7.22, 11.077, 5.625], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")),
        PharmacophoricPoint(
        feat_type="hb donor",
        center=puw.quantity([4.778, 8.432, 7.805], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")),
        PharmacophoricPoint(
        feat_type="aromatic ring",
        center=puw.quantity([1.56433333333334, 7.06399999999999, 3.135], "angstroms"),
        radius=puw.quantity(1.0, "angstroms"))
    ]
    pharmacophore = Pharmacophore(elements)

    file_path = "./openpharmacophore/data/ligands/mols.smi"

    screener = screening3D.VirtualScreening3D(pharmacophore)
    screener.screen_db_from_dir(file_path, titleLine=False)

    assert screener.n_molecules == 5
    assert screener.n_matches == 4
    assert screener.n_fails == 1
    assert len(screener.matches) == 4
    for _, id, mol in screener.matches:
        assert id is None
        assert isinstance(mol, Chem.Mol)

### Tests for VirtrualScreening2D class ###
def test_screen_db_from_dir_2D():
    file_path = "./openpharmacophore/data/ligands/mols.smi"
    mol = Chem.MolFromSmiles("Cc1cccc(c2n[nH]cc2c3ccc4ncccc4n3)n1")

    screener = screening2D.VirtualScreening2D(mol, similarity="tanimoto", sim_cutoff=0.6)
    screener.screen_db_from_dir(file_path, titleLine=False)
    assert screener.n_molecules == 5
    assert screener.n_matches == 1
    assert screener.n_fails == 4
    assert len(screener.matches) == 1
    assert screener.matches[0][0] == 1.0
    assert screener.matches[0][1] is None
    assert isinstance(screener.matches[0][2], Chem.Mol)
