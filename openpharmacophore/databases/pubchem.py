from openpharmacophore._private_tools.exceptions import FetchError, OpenPharmacophoreValueError
import pandas as pd
from tqdm.auto import tqdm
from io import StringIO
import json
import requests
import time
from typing import List, Dict, Tuple

base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

def _get_data(url: str, attempts: int = 5) -> requests.Response.content:
    """ Downloads data from a given url.
        
        Parameters
        ----------
        url : str
            url to fetch data from
        
        attempts : int
            number of times to try to download the data in case of failure

        Returns
        ----------
        requests.Response.content
            The content of the response
    """
    if attempts <= 0:
        raise ValueError("Number of attempts must be greater than 0")

    while (attempts > 0):
        res = requests.get(url)
        if res.status_code == requests.codes.ok:
            break
        else:
            attempts -= 1
            time.sleep(2)
    
    if res.status_code != requests.codes.ok:
        raise FetchError("Failed to get data from {}".format(url))
    
    return res.content
    

def get_assay_compounds_id(assay_id: int, attempts: int = 10) -> List[int]:
    """ Get compounds id for tested compounds in an assay

        Parameters
        ----------
        assay_id : int
            The id of the bioassay.
        
        attempts : int, optional
            Number of times to try to download the data in case of failure
            (default=10). 

        Returns
        -------
        list 
            List containing the compounds ids.
    """
    assay_url = base_url + "/bioassay/AID/{}/cids/JSON".format(assay_id)
    data = _get_data(assay_url, attempts)
    
    ids_dict = json.loads(data)

    return ids_dict["InformationList"]["Information"][0]["CID"]

def get_assay_description(assay_id: int, summary: bool =True, attempts: int = 10) -> Dict[str, str]:
    """ Get the description of an assay in JSON format.

        Parameters
        ----------
        assay_id : int
            The id of the bioassay.
        
        summary : bool
            If true returns a summary of the description of the assay (default=True).  

        attempts : int
            number of times to try to download the data in case of failure
            (default=10). 

        Returns
        --------
        dict
            A dictionary containing the assay description.
    """
    assay_url = base_url + "/assay/aid/{}".format(assay_id)

    if summary:
        description_url = assay_url + "/summary/JSON"
    else:
        description_url = assay_url + "/description/JSON"

    data = _get_data(description_url, attempts)
    return json.loads(data)

def get_assay_results(assay_id, form="dataframe", attempts=10):
    """ Get results of an assay. 

        Parameters
        ----------
        assay_id : int
            The id of the bioassay.
        
        form : {"dataframe", "dict"}
            The form of the returned object. Can be a dataframe or dict.  

        attempts : int, optional
            number of times to try to download the data in case of failure
            (default=10). 

        Returns
        ----------
        pandas.DataFrame or dict
            Dataframe or dictionarty with the assay results.
        
    """
    if form == "dataframe":
        format = "CSV"
    elif form == "dict":
        format = "JSON"
    else:
        raise ValueError("{} is not a valid form".format(form))

    assay_url = base_url + "/assay/aid/{}/{}".format(assay_id, format)

    data = _get_data(assay_url, attempts)
    
    if format == "CSV":
        csv_string = StringIO(data.decode("utf-8"))
        return pd.read_csv(csv_string)
    elif format == "JSON":
        return json.loads(data.content)
    
def get_assay_target_info(assay_id: int, attempts: int = 10) -> Dict[str, str]:
    """ Get target information of an assay.

        Parameters
        ----------
        assay_id : int
            The id of the bioassay.
        
        attempts : int, optional
            number of times to try to download the data in case of failure
            (default=10). 

        Returns
        --------
        dict
            A dictionary containing information of the assay target.
    """
    target_url = base_url + "/assay/aid/{}/targets/ProteinGI,ProteinName,GeneID,GeneSymbol/JSON".format(assay_id)
    data = _get_data(target_url, attempts)
    return json.loads(data)

def get_assay_bioactivity_data(assay_id):
    """ Get bioactivity data and the compounds in an assay. 
    
        Parameters
        ----------
        assay_id : int
            The id of the bioassay.
        Returns
        ----------
        compounds : list of 2-tuples
            The first element is the compound PubChem id.
            The second element is the smiles of the compound.
        
        bioactivity : np.array of bits
            An array where each element corresponds to the index of the compounds list.
            An entry is either one if the compound is active or zero if the compund is inactive.
    """
    assay_results = get_assay_results(assay_id=assay_id, form="dataframe")
    # Keep only cid and activity columns
    df = assay_results[["PUBCHEM_CID", "PUBCHEM_ACTIVITY_OUTCOME"]]
    df = df.dropna()
    # Add activity column of 1 and 0
    df["activity"] = df["PUBCHEM_ACTIVITY_OUTCOME"].apply(lambda x: 0 if x == "Inactive" else 1)
    # Drop original activity column
    df = df.drop("PUBCHEM_ACTIVITY_OUTCOME", axis=1)
    df = df.astype("int32")
    molecules_ids = df["PUBCHEM_CID"].tolist()
    bioactivity = df["activity"].to_numpy() 

    molecules = []
    print("Fetching molecules smiles...")
    for mol_id in tqdm(molecules_ids):
        smiles = get_compound_smiles(mol_id)
        molecules.append((mol_id, smiles))

    return molecules, bioactivity

def get_assay_actives_and_inactives(assay_id):
    """ Get smiles for compounds in an assay split into active and inactive.

        Parameters
        ----------
        assay_id : int
            The id of the bioassay.
        
        Returns
        ----------
        actives : 2-tuple
            The first element is a list of the active compounds PubChem ids, and
            the second elment is a list of smiles for the active compounds.
        
        inactives : 2-tuple
            The first element is a list of the inactive compounds PubChem ids, and
            the second elment is a list of smiles for the inactive compounds.
    """
    assay_results = get_assay_results(assay_id=assay_id, form="dataframe")
    # Keep only cid and activity columns
    df = assay_results[["PUBCHEM_CID", "PUBCHEM_ACTIVITY_OUTCOME"]]
    df = df.dropna()
    # Split into active/inactive
    actives = df.loc[df["PUBCHEM_ACTIVITY_OUTCOME"] == "Active"]
    inactives = df.loc[df["PUBCHEM_ACTIVITY_OUTCOME"] == "Inactive"]
    # Drop activity column
    actives = actives.drop("PUBCHEM_ACTIVITY_OUTCOME", axis=1)
    inactives = inactives.drop("PUBCHEM_ACTIVITY_OUTCOME", axis=1)
    # Cast to int
    actives = actives.astype("int32")
    inactives = inactives.astype("int32")

    actives_list = actives["PUBCHEM_CID"].tolist()
    inactives_list = inactives["PUBCHEM_CID"].tolist()
    
    actives_smiles = []
    inactives_smiles = []

    print("Fetching active compound smiles...")
    for compound in tqdm(actives_list):
        smiles = get_compound_smiles(compound)
        actives_smiles.append(smiles)
    
    print("Fetching inactive compound smiles...")
    for compound in tqdm(inactives_list):
        smiles = get_compound_smiles(compound)
        inactives_smiles.append(smiles)
            
    return (actives_list, actives_smiles), (inactives_list, inactives_smiles)

def get_compound_assay_summary(compound_id, form="dataframe", attempts=10):
    """ Get summary of biological test results for a given compound.

        Parameters
        ----------
        compound_id : int
            The PubChem id of the compound.
        
        form : {"dataframe", "dict"}
            The form of the returned object. Can be a dataframe or dict.

        attempts : int, optional
            Number of times to try to download the data in case of failure
            (default=10). 

        Returns
        ----------
        pandas.Dataframe or dict
            A dataFrame or a dictionary with the assay results for the passed compound.
    """
    if form == "dataframe":
        format = "CSV"
    elif form == "dict":
        format = "JSON"
    else:
        raise OpenPharmacophoreValueError("{} is not a valid form".format(form))
    
    compound_url = base_url + "/compound/cid/{}/assaysummary/{}".format(compound_id, format)
    data = _get_data(compound_url, attempts)
    
    if format == "CSV":
        csv_string = StringIO(data.decode("utf-8"))
        return pd.read_csv(csv_string)
    elif format == "JSON":
        return json.loads(data.content)

def get_compound_id(name: str, attempts: int = 10) -> str:
    """ Get pubchem compound id for a given compound name.

        Parameters
        ----------
        name : str
            Name of the compound. 
        
        attempts : int, optional
            Number of times to try to download the data in case of failure
            (default=10). 

        Returns
        ----------
        str
            A string with the compound name.
    """
    compound_url = base_url + "/compound/name/{}/cids/JSON".format(name)
    data = _get_data(compound_url, attempts)

    json_data = json.loads(data)
    return json_data["IdentifierList"]["CID"][0]

def get_compound_description(compound_identifier, attempts=10):
    """ Get description for a given compound. 

        Parameters
        ----------
        compound_identifier : str or int
            The name as str or the PubChem id as int of the compound.
        
        attempts : int, optional
            Number of times to try to download the data in case of failure
            (default=10). 

        Returns
        --------
        dict
            A dictionary containing the compound description.
    """
    # If a string is passed assume its compound name
    if isinstance(compound_identifier, str):
        compound_url = base_url + "/compound/name/{}/description/JSON".format(compound_identifier)
    # Else use compound id
    else:
        compound_url = base_url + "/compound/cid/{}/description/JSON".format(compound_identifier)

    data = _get_data(compound_url, attempts)
    return json.loads(data)

def get_compound_smiles(compound_id: int, attempts: int = 10) -> str:
    """ Get smiles for a given compound. 

        Parameters
        ----------
        compound_id: int
            The PubChem id of the compound.
        
        attempts : int, optional
            Number of times to try to download the data in case of failure
            (default=10). 

        Returns
        -------
        smiles : str 
            The smiles for the passed compound. 

    """
    smiles_url = base_url + "/compound/cid/{}/property/CanonicalSMILES/TXT".format(compound_id)
    data = _get_data(smiles_url, attempts)
    smiles = data.decode("utf-8").rstrip()
    return smiles

def get_target_assays(identifier: str, identifier_type: str, attempts: int = 10) -> pd.DataFrame:
    """ Get assay ids and name for a given target

        Parameters
        ----------
        identifier : str
            Identifier of the target. Can be GI, Gene ID, or gene symbol.

        identifer_type : {"genesymbol", "geneid", "gi"}
            The type of the identifier can be genesymbol, geneid or gi. 

        Returns
        --------
        pandas.DataFrame 
            Dataframe with the assays ids and names for the passed target
    """
    identifier_type = identifier_type.lower()
    valid_identifiers = ["genesymbol", "geneid", "gi"]
    if identifier_type not in valid_identifiers:
        raise OpenPharmacophoreValueError("{} is not a valid identifier type")

    target_url = base_url + "/assay/target/{}/{}/description/JSON".format(identifier_type, identifier)
    data = _get_data(target_url, attempts)

    assays_dict = json.loads(data)
    ids = []
    names = []
    # Retrieve only id and name
    for i in range(len(assays_dict["PC_AssayContainer"])):
        assay = assays_dict['PC_AssayContainer'][i]['assay']['descr']
        id = int(assay["aid"]["id"])
        name = assay["name"]
        ids.append(id)
        names.append(name)

    assays = {
        "id": ids,
        "name": names,
    }

    return pd.DataFrame.from_dict(assays)

def similarity_search(compound, threshold=None, max_records=None, attempts=5):
    """ Perform a 2D similarity search for a given compound.

        Parameters
        ----------
        compound : str or int
            Can be a smiles str or the PubChem id of the compound.
        
        threshold : int, optional
            Minimum Tanimoto score for a hit.
        
        max_records : int, optional
            Maximum number of hits.

        attempts : int, optional
            Number of times to try to download the data in case of failure
            (default=10). 

        Returns
        --------
        list
            A List with the id's of the hits.
    """
    # If compound is passed as str assume is smiles
    if isinstance(compound, str):
        url = base_url + "/compound/similarity/smiles/{}/JSON".format(compound)
    # Else use compound id
    else:
        url = base_url + "/compound/similarity/cid/{}/JSON".format(compound)
    
    if threshold and max_records:
        url += "?Threshold={}&MaxRecords={}".format(threshold, max_records)
    elif threshold:
        url += "?Threshold={}".format(threshold)
    elif max_records:
        url += "?MaxRecords={}".format(max_records)

    data = _get_data(url, attempts)
    # Data returns a listkey that can be used to retrieve the results from another url
    content = json.loads(data)
    listkey = content["Waiting"]["ListKey"]
    
    results_url = base_url + "/compound/listkey/{}/cids/JSON".format(listkey)
    # Wait a little as similarity searches take more time to complete
    time.sleep(5)
    data = _get_data(results_url, attempts)
    data_dict = json.loads(data)

    return data_dict["IdentifierList"]["CID"]

