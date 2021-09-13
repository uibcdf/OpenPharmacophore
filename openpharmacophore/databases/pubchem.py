from openpharmacophore._private_tools.exceptions import InvalidFileFormat
import pandas as pd
from tqdm.auto import tqdm
from io import StringIO
import json
import requests
import time

class PubChem():

    def __init__(self):
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    def _get_data(self, url, attempts):

        status_code = None
        while (attempts > 0):
            try:
                res = requests.get(url)
                status_code = int(res.status_code)
                if status_code == 200 or status_code == 202:
                    break
                else:
                    attempts -= 1
                    time.sleep(2)
            except Exception as e:
                print(e)
                attempts -= 1
                time.sleep(2)
        
        if status_code != 200 and status_code != 202:
            raise Exception("Failed to get data from {}".format(url))
        
        return res.content
        

    def get_assay_compounds_id(self, assay_id, attempts=10):
        """ Get compounds id for tested compunds in an assay

            Parameters
            ----------
            Returns
            ----------
        """
        assay_url = self.base_url + "/bioassay/AID/{}/cids/JSON".format(assay_id)
        data = self._get_data(assay_url, attempts)
        
        return json.loads(data)

    def get_assay_description(self, assay_id, summary=True, attempts=10):
        """ Get the description of an assay in JSON format

            Parameters
            ----------
            Returns
            ----------
        """
        assay_url = self.base_url + "/assay/aid/{}".format(assay_id)

        if summary:
            description_url = assay_url + "/summary/JSON"
        else:
            description_url = assay_url + "/description/JSON"

        data = self._get_data(description_url, attempts)
        return json.loads(data)

    def get_assay_results(self, assay_id, format="csv", attempts=10):
        """ Get results of an assay 

            Parameters
            ----------
            Returns
            ----------
        """
        format = format.upper()
        assay_url = self.base_url + "/assay/aid/{}/{}".format(assay_id, format)

        data = self._get_data(assay_url, attempts)
       
        if format == "CSV":
            csv_string = StringIO(data.decode("utf-8"))
            return pd.read_csv(csv_string)
        elif format == "JSON":
            return json.loads(data.content)
        
    def get_assay_target_info(self, assay_id, attempts=10):
        """ Get target information of an assay 

            Parameters
            ----------
            Returns
            ----------
        """
        target_url = self.base_url + "/assay/aid/{}/targets/ProteinGI,ProteinName,GeneID,GeneSymbol/JSON".format(assay_id)
        data = self._get_data(target_url, attempts)
        return json.loads(data)
    
    def get_assay_training_data(self, assay_id):
        """ Get smiles for compounds in an assay split into active and inactive

            Parameters
            ----------
            Returns
            ----------
        """
        assay_results = self.get_assay_results(assay_id=assay_id, format="csv")
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
            smiles = self.get_compound_smiles(compound)
            actives_smiles.append(smiles)
        
        print("Fetching inactive compound smiles...")
        for compound in tqdm(inactives_list):
            smiles = self.get_compound_smiles(compound)
            inactives_smiles.append(smiles)
                
        return (actives_list, actives_smiles), (inactives_list, inactives_smiles)

    def get_compound_assay_summary(self, compound_id, format="csv", attempts=10):
        """ Get summary of biological test results for a given compound

            Parameters
            ----------
            Returns
            ----------
        """
        format = format.upper()
        compound_url = self.base_url + "/compound/cid/{}/assaysummary/{}".format(compound_id, format)
        data = self._get_data(compound_url, attempts)
       
        if format == "CSV":
            csv_string = StringIO(data.decode("utf-8"))
            return pd.read_csv(csv_string)
        elif format == "JSON":
            return json.loads(data.content)
    
    def get_compound_id(self, name, attempts=10):
        """ Get pubchem compound id for a given compound name 

            Parameters
            ----------
            Returns
            ----------
        """
        compound_url = self.base_url + "/compound/name/{}/cids/JSON".format(name)
        data = self._get_data(compound_url, attempts)

        json_data = json.loads(data)
        return json_data["IdentifierList"]["CID"][0]

    def get_compound_description(self, compound_identifier, attempts=10):
        """ Get description for a given compound 

            Parameters
            ----------
            Returns
            ----------
        """
        # If a string is passed assume its compund name
        if isinstance(compound_identifier, str):
            compound_url = self.base_url + "/compound/name/{}/description/JSON".format(compound_identifier)
        # Else use compound id
        else:
            compound_url = self.base_url + "/compound/cid/{}/description/JSON".format(compound_identifier)

        data = self._get_data(compound_url, attempts)
        return json.loads(data)

    def get_compound_smiles(self, compound_id, attempts=10):
        """ Get smiles for a given compound 

            Parameters
            ----------
            Returns
            ----------
        """
        smiles_url = self.base_url + "/compound/cid/{}/property/CanonicalSMILES/TXT".format(compound_id)
        data = self._get_data(smiles_url, attempts)
        smiles = data.decode("utf-8").rstrip()
        return smiles

    def get_target_assays(self, identifier, identifier_type, attempts=10):
        """ Get assay ids for a given target

            Parameters
            ----------
            Returns
            ----------
        """
        identifier_type = identifier_type.lower()
        target_url = self.base_url + "assay/target/{}/{}/aids/JSON".format(identifier_type, identifier)
        data = self._get_data(target_url, attempts)
        return json.load(data)

    def similarity_search(self, compound, threshold=None, max_records=None, attempts=10):
        """ Perform a 2D similarity search for a given compound

            Parameters
            ----------
            Returns
            ----------
        """
        # If compound is passed as str assume is smiles
        if isinstance(compound, str):
            url = self.base_url + "/compound/similarity/smiles/{}/JSON".format(compound)
        # Else use compound id
        else:
            url = self.base_url + "/compound/similarity/cid/{}/JSON".format(compound)
        
        if threshold and max_records:
            url += "?Threshold={}&MaxRecords={}".format(threshold, max_records)
        elif threshold:
            url += "?Threshold={}".format(threshold)
        elif max_records:
            url += "?MaxRecords={}".format(max_records)

        data = self._get_data(url, attempts)
        # Data returns a listkey that can be used to retrieve the results from another url
        content = json.loads(data)
        listkey = content["Waiting"]["ListKey"]
        
        results_url = self.base_url + "/compound/listkey/{}/cids/JSON".format(listkey)
        # Wait a little as similarity searches take more time to complete
        time.sleep(5)
        data = self._get_data(results_url, attempts)

        return json.loads(data)

