import openpharmacophore as oph
import openpharmacophore.data as data
from screening_pharmacophores import four_element_pharmacophore, pharmacophore_fingerprint
import numpy as np
import pandas as pd
import pytest
from collections import namedtuple


def test_retrospective_screening_constructor():
    pharmacophore = four_element_pharmacophore()
    fingerprint = pharmacophore_fingerprint()

    screener_3d = oph.RetrospectiveScreening(pharmacophore)
    assert screener_3d.scoring_metric == "SSD"
    assert screener_3d._screen_fn.__name__ == "_align_molecules"

    screener_2d = oph.RetrospectiveScreening(fingerprint, similarity="dice")
    assert screener_2d.scoring_metric == "Similarity"
    assert screener_2d._screen_fn.__name__ == "_fingerprint_similarity"


def test_from_bioactivity_data():

    molecules = [
        (663426, 'CCOC1=CC=CC(=C1)C2=NN=C3N2N=C(S3)C4=CC(=C(C(=C4)OC)OC)OC'),
        (9550244, 'CC1=CC(=NN1CC2=C(C=CC(=C2)C=C3CCC4=C(C3=O)C=CC(=C4)OC)OC)C'),
        (5770444, 'CC1=CC=C(C=C1)C2=NN=C(O2)CSC3=NC4=C(C5=C(N4)C=CC(=C5)OC)N=N3'),
        (644597, 'C1CCN(CC1)S(=O)(=O)C2=CC=C(C=C2)C3=NC(=NN3)SCC(=O)NC4=NC=CS4'),
        (644708, 'C1=CC=C2C(=C1)C=C(C=N2)NC(=O)C3=CC=CO3')
    ]
    bioactivity = np.array([0, 0, 1, 0, 0])
    pharmacophore = oph.StructuredBasedPharmacophore.from_pdb(
        data.pdb["1qku"], ligand_id="EST:A:600")
    pharmacophore.remove_points([0, 4])
    screener = oph.RetrospectiveScreening(pharmacophore)
    screener.from_bioactivity_data(molecules, bioactivity)
    assert screener.n_molecules == 5
    assert screener.n_actives == 1
    assert screener.n_inactives == 4

    assert screener.molecules[0].score == 0.0
    assert screener.molecules[1].score == 0.0
    assert round(screener.molecules[2].score, 2) == 0.08
    assert round(screener.molecules[3].score, 2) == 0.12
    assert screener.molecules[4].score == 0.0


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
                         1, 0, 0, 0, 0, ])

    scores_and_labels_2 = (scores_2, labels_2)

    return scores_and_labels_1, scores_and_labels_2


def test_auc(scores_and_labels):
    scores_and_labels_1, scores_and_labels_2 = scores_and_labels

    scores_1, labels_1 = scores_and_labels_1
    auc = oph.RetrospectiveScreening._get_auc(scores_1, labels_1)
    assert auc == 0.68

    scores_2, labels_2 = scores_and_labels_2
    auc = oph.RetrospectiveScreening._get_auc(scores_2, labels_2)
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

    fpr, tpr = oph.RetrospectiveScreening._roc_points(scores_1, labels_1)
    assert fpr == expected_fpr
    assert tpr == expected_tpr

    scores_2, labels_2 = scores_and_labels_2
    expected_fpr = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0]
    expected_tpr = [0.0, 0.3333333333333333, 0.5, 0.6666666666666666,
                    0.8333333333333334, 1.0, 1.0, 1.0, 1.0, 1.0]

    fpr, tpr = oph.RetrospectiveScreening._roc_points(scores_2, labels_2)
    assert fpr == expected_fpr
    assert tpr == expected_tpr


@pytest.fixture
def dataset_for_enrichment_tests():
    """ Load a dataset of the egfr bioassay with the bioactivity of each
        molecule as well as the scores of maccs and morgan fingerprints.
    """
    dataset = pd.read_csv(data.bioassays["egfr_bioassay"])
    maccs = dataset["tanimoto_maccs"].to_numpy()
    morgan = dataset["tanimoto_morgan"].to_numpy()
    bioactivity = dataset["activity"].to_numpy()

    return maccs, morgan, bioactivity


def test_enrichment_data(dataset_for_enrichment_tests):
    maccs_score, morgan_score, bioactivity = dataset_for_enrichment_tests

    # Tests for enrichment data with maccs fingerprints
    enrichment_df = pd.read_csv(data.bioassays["enrichment_data"])
    screen_percent_expected = enrichment_df["macss_per_screen"].to_numpy()
    screen_percent, actives_percent = oph.RetrospectiveScreening._enrichment_data(maccs_score, bioactivity)
    screen_percent = np.array(screen_percent)
    assert len(screen_percent) == len(screen_percent_expected)
    assert np.allclose(screen_percent, screen_percent_expected)
    assert screen_percent.shape[0] == len(actives_percent)

    # Tests for enrichment data with morgan fingerprints
    screen_percent_expected = enrichment_df["morgan_per_screen"].to_numpy()
    screen_percent, actives_percent = oph.RetrospectiveScreening._enrichment_data(morgan_score, bioactivity)
    screen_percent = np.array(screen_percent)
    assert len(screen_percent) == len(screen_percent_expected)
    assert np.allclose(screen_percent, screen_percent_expected)
    assert screen_percent.shape[0] == len(actives_percent)


def test_calculate_enrichment_factor(dataset_for_enrichment_tests):
    maccs_score, morgan_score, bioactivity = dataset_for_enrichment_tests

    # Maccs enrichment factor
    expected_enrichment = 7.318982387475538
    enrichment_factor = oph.RetrospectiveScreening._calculate_enrichment_factor(maccs_score, bioactivity, 5)
    assert expected_enrichment == enrichment_factor
    # Morgan enrichment factor
    expected_enrichment = 7.9060665362035225
    enrichment_factor = oph.RetrospectiveScreening._calculate_enrichment_factor(morgan_score, bioactivity, 5)
    assert expected_enrichment == enrichment_factor


def test_ideal_enrichment_factor(dataset_for_enrichment_tests):
    # Mock pharmacophore
    pharmacophore = pharmacophore_fingerprint()
    vs = oph.RetrospectiveScreening(pharmacophore)

    _, _, bioactivity = dataset_for_enrichment_tests
    vs.bioactivities = bioactivity
    vs.n_actives = np.sum(bioactivity)

    ideal_ef = vs.ideal_enrichment_factor(5)
    assert ideal_ef == 8.855185909980431


def test_confusion_matrix():

    screener = oph.RetrospectiveScreening(pharmacophore_fingerprint())
    screener.bioactivities = np.array([0, 1, 0, 1, 1, 0, 1])

    MolScoreMock = namedtuple("MolScore", ["score", "id", "mol"])
    molecules = [
        MolScoreMock(0.5, 1, ""),
        MolScoreMock(0.9, 1, ""),
        MolScoreMock(0.7, 1, ""),
        MolScoreMock(0.7, 1, ""),
        MolScoreMock(0.3, 1, ""),
        MolScoreMock(0.4, 1, ""),
        MolScoreMock(0.5, 1, ""),
    ]
    screener.molecules = molecules
    screener.n_molecules = len(molecules)

    threshold = 0.6

    expected_matrix = np.array(
        [[2, 1],
         [2, 2]]
    )

    assert np.all(screener.confusion_matrix(threshold) == expected_matrix)
