from openpharmacophore import PharmacophoricPoint, Pharmacophore, VirtualScreening, RetrospectiveScreening, pharmacophore
import numpy as np
import pytest
import pyunitwizard as puw
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D
from rdkit.Chem.Pharm2D.Generate import Gen2DFingerprint

### Tests for VitualScreening class with standard pharmacophore ###

@pytest.fixture
def mock_screening_results():
    """ Returns an isntance of a VirtualScreening class
        with fake attributes.
    """
    pharmacophore = Pharmacophore()
    screener = VirtualScreening(pharmacophore)
    screener.n_molecules = 60000
    screener.n_matches = 3
    screener.n_fails = 60000 - 3
    screener.scoring_metric = "Similarity"

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

           ID        Similarity
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
        assert "Id" in results
        assert results["Id"] == ["CHEMBL22796",
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
        assert list(results.columns) == ['Id', 'Smiles', 'Similarity', 'Mol_weight', 'logP']
        assert results["Id"].to_list() == ["CHEMBL22796",  "CHEMBL431887", "CHEMBL907"]
        assert results["Smiles"].to_list() == ["CCCCc1nn(CCC)c(C(=O)O)c1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1",
            "CCCCc1nn(CC)c(C(=O)O)c1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1",
            "CCCCc1nc(Cl)c(C(=O)O)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1"
        ]
        assert np.allclose(results["Mol_weight"].to_numpy(), np.array([444.54, 430.51, 436.90]))
        assert np.allclose(results["logP"].to_numpy(), np.array([4.7718, 4.3817, 4.4727]))
        assert np.allclose(results["Similarity"].to_numpy(), np.array([0.6542, 0.7833, 0.8974]))

@pytest.fixture
def four_element_pharmacophore():
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
    return Pharmacophore(elements)


@pytest.fixture
def pharmacophore_fingerprint():
    mol = Chem.MolFromSmiles("Cc1cccc(c2n[nH]cc2c3ccc4ncccc4n3)n1")
    factory = Gobbi_Pharm2D.factory
    return Gen2DFingerprint(mol, factory)


def test_screen_mol_file(pharmacophore_fingerprint, four_element_pharmacophore):
    file_path = "./openpharmacophore/data/ligands/mols.smi"
   
    screener = VirtualScreening(pharmacophore_fingerprint, similarity="tanimoto", sim_cutoff=0.6)
    screener.screen_mol_file(file_path)
    assert screener.n_molecules == 5
    assert screener.n_matches == 1
    assert screener.n_fails == 4
    assert len(screener.matches) == 1
    assert screener.matches[0][0] == 1.0
    assert screener.matches[0][1] is None
    assert isinstance(screener.matches[0][2], Chem.Mol)
    
    pharmacophore = four_element_pharmacophore
    screener = VirtualScreening(pharmacophore)
    screener.screen_mol_file(file_path)

    assert screener.n_molecules == 5
    assert screener.n_matches == 4
    assert screener.n_fails == 1
    assert len(screener.matches) == 4
    for _, id, mol in screener.matches:
        assert id is None
        assert isinstance(mol, Chem.Mol)



## Tests for Retrospective Screening ##
def test_init_retrospective_screening(four_element_pharmacophore, pharmacophore_fingerprint):
    
    pharmacophore = four_element_pharmacophore
    fingerprint = pharmacophore_fingerprint
    
    screener_3d = RetrospectiveScreening(pharmacophore)
    assert screener_3d.scoring_metric == "SSD"
    assert screener_3d._screen_fn.__name__ == "_align_molecules"
    
    screener_2d = RetrospectiveScreening(fingerprint, similarity="dice")
    assert screener_2d.scoring_metric == "Similarity"
    assert screener_2d._screen_fn.__name__ == "_fingerprint_similarity"

def test_from_bioactivity_data():
    pass


@pytest.fixture
def scores_and_labels():
    # Example test set 1  
    scores_1 = np.array([0.90, 0.80, 0.70, 0.60, 0.55,
                      0.54, 0.53, 0.52, 0.51, 0.505,
                      0.40, 0.39, 0.38, 0.37, 0.36,
                      0.35, 0.34, 0.33, 0.30, 0.10])

    labels_1 = np.array([1, 1, 0, 1, 1,
                       1, 0, 0, 1, 0,
                       1, 0, 1, 0, 0,
                       0, 1, 0, 1, 0])
    
    scores_and_labels_1 = (scores_1, labels_1)
    
    scores_2 = np.array([0.99999, 0.99999, 0.99993, 0.99986, 0.99964,
                       0.99955, 0.68139, 0.50961, 0.48880, 0.44951])

    labels_2 = np.array([1, 1, 1, 1, 1,
                       1, 0, 0, 0, 0,])
    
    scores_and_labels_2 = (scores_2, labels_2)
    
    return scores_and_labels_1, scores_and_labels_2

def test_auc(scores_and_labels):
    scores_and_labels_1, scores_and_labels_2 = scores_and_labels
    
    scores_1, labels_1 = scores_and_labels_1
    auc = RetrospectiveScreening._get_auc(scores_1, labels_1)
    assert auc == 0.68
    
    scores_2, labels_2 = scores_and_labels_2
    auc = RetrospectiveScreening._get_auc(scores_2, labels_2)
    assert auc == 1.0

def test_roc_points(scores_and_labels):
    scores_and_labels_1, scores_and_labels_2 = scores_and_labels
    
    scores_1, labels_1 = scores_and_labels_1
    expected_fpr = [0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.1, 
                    0.2, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5, 
                    0.6, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0]
    expected_tpr = [0.0, 0.1, 0.2, 0.2, 0.3, 0.4, 0.5, 
                    0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 
                    0.8, 0.8, 0.8, 0.9, 0.9, 1.0, 1.0]
    
    fpr, tpr = RetrospectiveScreening._roc_points(scores_1, labels_1)
    assert fpr == expected_fpr
    assert tpr == expected_tpr
    
    scores_2, labels_2 = scores_and_labels_2
    expected_fpr = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0]
    expected_tpr = [0.0, 0.3333333333333333, 0.5, 0.6666666666666666, 
                    0.8333333333333334, 1.0, 1.0, 1.0, 1.0, 1.0]
    
    fpr, tpr = RetrospectiveScreening._roc_points(scores_2, labels_2)
    assert fpr == expected_fpr
    assert tpr == expected_tpr

@pytest.fixture
def dataset_for_enrichment_tests():
    """ Load a dataset of the egfr bioassay with the bioactivity of each
        molecule as well as the scores of maccs and morgan fingerprints.
    """
    dataset = pd.read_csv("./openpharmacophore/data/bioassays/egfr_bioassay.csv")
    maccs = dataset["tanimoto_maccs"].to_numpy()
    morgan = dataset["tanimoto_morgan"].to_numpy()
    bioactivity = dataset["activity"].to_numpy()
    
    return maccs, morgan, bioactivity
    
def test_enrichment_data(dataset_for_enrichment_tests):
    maccs_score, morgan_score, bioactivity = dataset_for_enrichment_tests
    
    # Tests for enrichment data with maccs fingerprints
    enrichment_df = pd.read_csv("./openpharmacophore/data/bioassays/enrichment_data.csv")
    screen_percent_expected = enrichment_df["macss_per_screen"].to_numpy()
    screen_percent, actives_percent = RetrospectiveScreening._enrichment_data(maccs_score, bioactivity)
    screen_percent = np.array(screen_percent)
    assert len(screen_percent) == len(screen_percent_expected)
    assert np.allclose(screen_percent, screen_percent_expected)
    assert screen_percent.shape[0] == len(actives_percent)
    
    # Tests for enrichment data with morgan fingerprints
    screen_percent_expected = enrichment_df["morgan_per_screen"].to_numpy()
    screen_percent, actives_percent = RetrospectiveScreening._enrichment_data(morgan_score, bioactivity)
    screen_percent = np.array(screen_percent)
    assert len(screen_percent) == len(screen_percent_expected)
    assert np.allclose(screen_percent, screen_percent_expected)
    assert screen_percent.shape[0] == len(actives_percent)
    
def test__calculate_enrichment_factor(dataset_for_enrichment_tests):
    maccs_score, morgan_score, bioactivity = dataset_for_enrichment_tests
    
    # Maccs enrichment factor
    expected_enrichment = 7.318982387475538
    enrichment_factor = RetrospectiveScreening._calculate_enrichment_factor(maccs_score, bioactivity, 5)
    assert expected_enrichment == enrichment_factor
    # Morgan enrichment factor
    expected_enrichment = 7.9060665362035225
    enrichment_factor = RetrospectiveScreening._calculate_enrichment_factor(morgan_score, bioactivity, 5)
    assert expected_enrichment == enrichment_factor
    

def test_ideal_enrichment_factor(pharmacophore_fingerprint, dataset_for_enrichment_tests):
    # Mock pharmacophore
    pharmacophore = pharmacophore_fingerprint
    vs = RetrospectiveScreening(pharmacophore)
    
    _, _, bioactivity = dataset_for_enrichment_tests
    vs.bioactivities = bioactivity
    vs.n_actives = np.sum(bioactivity)
    
    ideal_ef = vs.ideal_enrichment_factor(5)
    assert ideal_ef == 8.855185909980431

def test_confussion_matrix():
    pass

