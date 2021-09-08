from openpharmacophore.databases.zinc import download_ZINC2D_smiles, discretize_values
import pytest

@pytest.mark.parametrize("subset,mol_weight,logp", [
    ("Drug-Like", None, None),
    (None, (250, 350), (-1, 1)),
    (None, (365, 415), (1.5, 2.25))
])
def test_download_ZINC2D_smiles(subset, mol_weight, logp):

    base_url = "http://files.docking.org/2D/"
    url_list = download_ZINC2D_smiles("./", 
        subset=subset,
        mol_weight_range=mol_weight, 
        logp_range=logp,
        testing=True)

    if subset == "Drug-like":
        assert len(url_list) == 90 * 4 * 2 
        assert url_list[0] == base_url + "BA/BAAA.smi"
        assert url_list[-1] == base_url + "JJ/JJEB.smi"
    elif mol_weight == (250, 350):
        assert len(url_list) == 12 * 4 * 2
        assert url_list[0] == base_url + "BA/BAAA.smi"
        assert url_list[-1] == base_url + "EC/ECEB.smi"
    elif mol_weight == (365, 415):
        assert len(url_list) == 12 * 4 * 2
        assert url_list[0] == base_url + "EC/ECAA.smi"
        assert url_list[-1] == base_url + "HE/HEEB.smi"

@pytest.mark.parametrize("value,lower", [
    (230, True),
    (484, False),
    (600, True)
])
def test_discretize_values(value, lower):
    
    bins = [200, 250, 300, 325, 350, 375, 400, 425, 450, 500, 550]
    new_value = discretize_values(value=value, bins=bins, name="Test", lower=lower)

    if value == 230:
        assert new_value == 200
    elif value == 484:
        assert new_value == 500
    else:
        assert new_value == 550