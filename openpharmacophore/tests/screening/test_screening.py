import openpharmacophore as oph
import openpharmacophore.data as data
from screening_pharmacophores import four_element_pharmacophore, pharmacophore_fingerprint
import numpy as np
import pytest
import pandas as pd
from rdkit import Chem


@pytest.fixture
def mock_screening_results():
    """ Returns an instance of a VirtualScreening class
        with fake attributes.
    """
    pharmacophore = oph.Pharmacophore()
    screener = oph.VirtualScreening(pharmacophore)
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
def test_get_screening_results(form, mock_screening_results):
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
        assert results["Id"].to_list() == ["CHEMBL22796", "CHEMBL431887", "CHEMBL907"]
        assert results["Smiles"].to_list() == ["CCCCc1nn(CCC)c(C(=O)O)c1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1",
                                               "CCCCc1nn(CC)c(C(=O)O)c1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1",
                                               "CCCCc1nc(Cl)c(C(=O)O)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1"
                                               ]
        assert np.allclose(results["Mol_weight"].to_numpy(), np.array([444.54, 430.51, 436.90]))
        assert np.allclose(results["logP"].to_numpy(), np.array([4.7718, 4.3817, 4.4727]))
        assert np.allclose(results["Similarity"].to_numpy(), np.array([0.6542, 0.7833, 0.8974]))


def test_screen_mol_file():
    file_path = data.ligands["mols"]

    screener = oph.VirtualScreening(pharmacophore_fingerprint(), similarity="tanimoto", sim_cutoff=0.6)
    screener.screen_mol_file(file_path)
    assert screener.n_molecules == 5
    assert screener.n_matches == 1
    assert screener.n_fails == 4
    assert len(screener.matches) == 1
    assert screener.matches[0][0] == 1.0
    assert screener.matches[0][1] is None
    assert isinstance(screener.matches[0][2], Chem.Mol)

    pharmacophore = four_element_pharmacophore()
    screener = oph.VirtualScreening(pharmacophore)
    screener.screen_mol_file(file_path)

    assert screener.n_molecules == 5
    assert screener.n_matches == 4
    assert screener.n_fails == 1
    assert len(screener.matches) == 4
    for _, id, mol in screener.matches:
        assert id is None
        assert isinstance(mol, Chem.Mol)
