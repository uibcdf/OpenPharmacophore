from openpharmacophore.databases.zinc_client import ZincClient
from openpharmacophore._private_tools import exceptions as exc
import pytest

def test_validate_filters():
    
    client = ZincClient()
    
    with pytest.raises(exc.InvalidFileFormat, match="is not a valid fileformat"):
        client._validate_filters(fileformat="pdf")
    
    with pytest.raises(exc.InvalidAvailabilityError, match="is not a valid availability"):
        client._validate_filters(fileformat="smi", availability="instantaneous")
        
    with pytest.raises(exc.InvalidBioactiveError, match="is not a valid bioactivity"):
        client._validate_filters(fileformat="smi", bioactive="hello")
    
    with pytest.raises(exc.InvalidBiogenicError, match="is not a valid biogenic"):
        client._validate_filters(fileformat="smi", biogenic="hello")
        
    with pytest.raises(exc.InvalidReactivityError, match="is not a valid reactivity"):
        client._validate_filters(fileformat="smi", reactivity="explosive")
    
def test_append_filters_to_url():
    
    client = ZincClient()
    catalog_url = client._catalog_url + "/chembl27/substances"
    
    full_url = client._append_filters_to_url(catalog_url, "smi", count=1000)
    assert full_url == "https://zinc.docking.org/catalogs/chembl27/substances.smi?count=1000"
    
    full_url = client._append_filters_to_url(catalog_url, "smi", count=1000, availability="for-sale")
    assert full_url == "https://zinc.docking.org/catalogs/chembl27/substances/subsets/for-sale.smi?count=1000"
    
    full_url = client._append_filters_to_url(catalog_url, "smi", count=1000, availability="for-sale", reactivity="anodyne")
    assert full_url == "https://zinc.docking.org/catalogs/chembl27/substances/subsets/for-sale+anodyne.smi?count=1000"
    
    full_url = client._append_filters_to_url(catalog_url, "smi", count=1000, bioactive="fda")
    assert full_url == "https://zinc.docking.org/catalogs/chembl27/substances/subsets/fda.smi?count=1000"
    
    full_url = client._append_filters_to_url(catalog_url, "smi", count=1000, biogenic="metabolites")
    assert full_url == "https://zinc.docking.org/catalogs/chembl27/substances/subsets/metabolites.smi?count=1000"
    
    substances_url = client._substances_url
    
    full_url = client._append_filters_to_url(substances_url, "smi", count=1000)
    assert full_url == "https://zinc.docking.org/substances.smi?count=1000"
    
    full_url = client._append_filters_to_url(substances_url, "smi", count=1000, availability="for-sale", reactivity="anodyne")
    assert full_url == "https://zinc.docking.org/substances/subsets/for-sale+anodyne.smi?count=1000"
    

def test_get_catalog_url():
    
    client = ZincClient()
    file_name = "./temp.smi"
    
    # First we test for errors
    catalog_name = "hello"
    with pytest.raises(exc.InvalidCatalogError, match="is not a valid catalog name"):
        client._get_catalog_url(file_name, catalog_name)
    
    catalog_name = "Asinex"
    assert client._get_catalog_url(file_name, catalog_name, count="all") == "https://zinc.docking.org/catalogs/asin/substances.smi?count=all"
    
    catalog_name = "Apollo Scientific"
    catalog_url = client._get_catalog_url(file_name, catalog_name, availability="for-sale") 
    assert catalog_url == "https://zinc.docking.org/catalogs/apollo/substances/subsets/for-sale.smi?count=1000"
    
    catalog_name = "ChemBridge Economical"
    catalog_url = client._get_catalog_url(file_name, catalog_name, availability="wait-ok", reactivity="clean") 
    assert catalog_url == "https://zinc.docking.org/catalogs/chbre/substances/subsets/wait-ok+clean.smi?count=1000"
    

def test_urls_for_tranches_2d():
    
    client = ZincClient()
    col_list = ["A", "B", "C"]
    row_list = ["A", "B", "C"]
    
    url_list = client._urls_for_tranches_2d(col_list, row_list)
    assert len(url_list) == 144
    # Test the first four
    assert url_list[0] == "http://files.docking.org/2D/AA/AAAA.smi"
    assert url_list[1] == "http://files.docking.org/2D/AA/AAAB.smi"
    assert url_list[2] == "http://files.docking.org/2D/AA/AAAC.smi"
    assert url_list[3] == "http://files.docking.org/2D/AA/AAAD.smi"
    # Test the last four
    assert url_list[-1] == "http://files.docking.org/2D/CC/CCED.smi"
    assert url_list[-2] == "http://files.docking.org/2D/CC/CCEC.smi"
    assert url_list[-3] == "http://files.docking.org/2D/CC/CCEB.smi"
    assert url_list[-4] == "http://files.docking.org/2D/CC/CCEA.smi"

    col_list = ["A", "B", "C", "D"]
    row_list = ["A", "B", "C", "D"]
    
    url_list = client._urls_for_tranches_2d(col_list, row_list)
    assert len(url_list) == 256 # len(col_list) * len(row_list) * 16
    
    col_list = ["B", "C"]
    row_list = ["A", "B", "C", "D"]
    
    url_list = client._urls_for_tranches_2d(col_list, row_list)
    assert len(url_list) == 128 # len(col_list) * len(row_list) * 16

def test_predifined_subset_tranches():
    
    client = ZincClient()
    
    # Test for errors
    with pytest.raises(exc.InvalidSubsetError, match="is not a valid subset"):
        client._predefined_subset_tranches(subset="hello")
    
    # Test each subset
    subset="Drug-Like"
    col_list, row_list = client._predefined_subset_tranches(subset)
    assert col_list == ["B", "C", "D", "E", "F", "G", "H", "I", "J"]
    assert row_list == ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
    
    subset="Lead-Like"
    col_list, row_list = client._predefined_subset_tranches(subset)
    assert col_list == ["C", "D", "E"]
    assert row_list == ["A", "B", "C", "D", "E", "F", "G"]
    
    subset="Lugs"
    col_list, row_list = client._predefined_subset_tranches(subset)
    assert col_list == ["E", "F", "G", "H", "I"]
    assert row_list == ["A", "B", "C", "D", "E", "F", "G", "H", "I"]
    
    subset="Goldilocks"
    col_list, row_list = client._predefined_subset_tranches(subset)
    assert col_list == ["C", "D", "E"]
    assert row_list == ["D", "E", "F"]
    
    subset="Big-n-Greasy"
    col_list, row_list = client._predefined_subset_tranches(subset)
    assert col_list == ["J", "K"]
    assert row_list == ["I", "J", "K"]

    subset="Fragments"
    col_list, row_list = client._predefined_subset_tranches(subset)
    assert col_list == ["A", "B"]
    assert row_list == ["A", "B", "C", "D", "E", "F", "G"]
    
    subset="Flagments"
    col_list, row_list = client._predefined_subset_tranches(subset)
    assert col_list == ["B", "C", "D"]
    assert row_list == ["A", "B", "C", "D", "E", "F", "G"]
    
    subset="Shards"
    col_list, row_list = client._predefined_subset_tranches(subset)
    assert col_list == ["A"]
    assert row_list == ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"]

def test_molw_and_logp_tranches():
    
    client = ZincClient()
    mw = (250, 200)
    logp = (0, 2)
    
    #Test for errors first
    with pytest.raises(exc.InvalidMolecularWeightRangeError, match="First number must be smaller"):
        client._mw_and_logp_tranches(mw, logp)
    
    mw = (100, 200)
    with pytest.raises(exc.InvalidMolecularWeightRangeError, match="Molecular weight must be a value between"):
        client._mw_and_logp_tranches(mw, logp)
        
    mw = (300, 700)
    with pytest.raises(exc.InvalidMolecularWeightRangeError, match="Molecular weight must be a value between"):
        client._mw_and_logp_tranches(mw, logp)
    
    mw = (250, 325)
    logp = (3, 1)
    with pytest.raises(exc.InvalidLogPRangeError, match="First number must be smaller"):
        client._mw_and_logp_tranches(mw, logp)
        
    logp = (-5, 2)
    with pytest.raises(exc.InvalidLogPRangeError, match="LogP must be a value between"):
        client._mw_and_logp_tranches(mw, logp)

    logp = (3, 10)
    with pytest.raises(exc.InvalidLogPRangeError, match="LogP must be a value between"):
        client._mw_and_logp_tranches(mw, logp)

    # Test that it returns the correct tranches
    mw = (200, 550)
    logp = (-1, 6)
    col_list, row_list = client._mw_and_logp_tranches(mw, logp)
    assert col_list == ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"]
    assert row_list == ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"]
    
    mw = (320, 500)
    logp = (0.5, 3)
    col_list, row_list = client._mw_and_logp_tranches(mw, logp)
    assert col_list == ["C", "D", "E", "F", "G", "H", "I", "J"]
    assert row_list == ["B", "C", "D", "E", "F"]
    
    mw = (415.5, 525)
    logp = (2.75, 4)
    col_list, row_list = client._mw_and_logp_tranches(mw, logp)
    assert col_list == ["G", "H", "I", "J", "K"]
    assert row_list == ["E", "F", "G", "H"]
    

def test_tranch_url_with_filters():
    
    client = ZincClient()
    url = client._tranch_url_with_filters(tranch="AAAA",
                                          availability=None,
                                          bioactive=None,
                                          biogenic=None,
                                          reactivity=None)
    
    assert url == 'https://zinc.docking.org/substances.smi?count=all&tranche_name=AAAA'
    
    url = client._tranch_url_with_filters(tranch="BBAA",
                                          availability="for-sale",
                                          bioactive=None,
                                          biogenic=None,
                                          reactivity="anodyne")
    
    assert url == "https://zinc.docking.org/substances/subsets/for-sale+anodyne.smi?count=all&tranche_name=BBAA"
    
    url = client._tranch_url_with_filters(tranch="CAAB",
                                          availability="for-sale",
                                          bioactive=None,
                                          biogenic="biogenic",
                                          reactivity=None,
                                          fileformat="sdf")
    
    assert url == "https://zinc.docking.org/substances/subsets/for-sale+biogenic.sdf?count=all&tranche_name=CAAB"
    
    url = client._tranch_url_with_filters(tranch="JCAA",
                                          availability="now",
                                          bioactive="in-vitro",
                                          biogenic=None,
                                          reactivity=None,
                                          fileformat="json"
                                          )

    assert url == "https://zinc.docking.org/substances/subsets/now+in-vitro.json?count=all&tranche_name=JCAA"

def test_tranche_with_filters_url_list():
    
    client = ZincClient()
    col_list = ["C"]
    row_list = ["B"]
    
    url_list = client._tranche_with_filters_url_list(col_list, row_list,
                                                     availability="for-sale",
                                                     bioactive=None,
                                                     biogenic=None,
                                                     reactivity=None,
                                                     fileformat="smi")
    assert len(url_list) == 16
    assert url_list[0] == "https://zinc.docking.org/substances/subsets/for-sale.smi?count=all&tranche_name=CBAA"
    assert url_list[-1] == "https://zinc.docking.org/substances/subsets/for-sale.smi?count=all&tranche_name=CBED"
    
    col_list = ["E", "F"]
    row_list = ["B"]
    
    url_list = client._tranche_with_filters_url_list(col_list, row_list,
                                                     availability="now",
                                                     bioactive=None,
                                                     biogenic=None,
                                                     reactivity="anodyne",
                                                     fileformat="mol2")
    assert len(url_list) == 32
    assert url_list[0] == "https://zinc.docking.org/substances/subsets/now+anodyne.mol2?count=all&tranche_name=EBAA"
    assert url_list[-1] == "https://zinc.docking.org/substances/subsets/now+anodyne.mol2?count=all&tranche_name=FBED"
    

def test_predifined_subset_3d_urls():
    pass

def test_molw_and_logp_3d_urls():
    pass

@pytest.mark.parametrize("value,lower", [
    (230, True),
    (484, False),
    (600, True)
])
def test_discretize_values(value, lower):
    
    bins = [200, 250, 300, 325, 350, 375, 400, 425, 450, 500, 550]
    new_value = ZincClient.discretize_values(value=value, bins=bins, name="Test", lower=lower)

    if value == 230:
        assert new_value == 200
    elif value == 484:
        assert new_value == 500
    else:
        assert new_value == 550
        