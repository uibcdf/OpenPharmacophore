from openpharmacophore._private_tools.exceptions import (InvalidLogPRangeError, InvalidZincIdError, InvalidCatalogError, InvalidFileFormat,
                        InvalidBioactiveError, InvalidBiogenicError, InvalidReactivityError,
                        InvalidAvailabilityError, CountTypeError, NegativeCountError,
                        InvalidSubsetError, OpenPharmacophoreTypeError, OpenPharmacophoreValueError, InvalidMolecularWeightRangeError, DownloadError,
                        ZincNotFoundError, ZincTimeoutError, InvaludUrlTypeError)
# Third Party
import requests
from tqdm.auto import tqdm
# Standard Library
import itertools
import json
import os
import pkg_resources
import re
import string

class ZincClient():
    """Class to interact with ZINC database."""
    
    def __init__(self):
        self._base_url = "https://zinc.docking.org/"
        self._substances_url = self._base_url + "substances"
        self._catalog_url = self._base_url + "catalogs"
        self._tranches_2d_url = "http://files.docking.org/2D/"

        self.catalogs = self._load_catalogs_from_file()            
        self.file_formats = ["xml" ,"csv","js",
                             "json","txt","mol2",
                             "db"  ,"sdf","smi",
                             "solv","db2"]
        
        self.filters = self._get_filters()
    
    ### Methods for downloading ###
    
    def compound_smiles(self, zinc_id):
        """ Get the smiles for a compound from its ZINC id.
        
            Parameters
            ----------
            zinc_id : str
                The zinc id of the compound.
                
            Returns
            -------
            str
                The smiles of the compound
        """
        if not isinstance(zinc_id, str):
            raise InvalidZincIdError(f"Invalid type {type(zinc_id)} Expected str.")
        
        pattern = r"ZINC[0-9]+"
        if re.match(pattern, zinc_id) is None:
            raise InvalidZincIdError(f"{zinc_id} is not a valid ZINC id")
        
        url = self._substances_url + "/" + zinc_id + ".json"
        res = requests.get(url, allow_redirects=True)
        if res.status_code != requests.codes.ok:
            raise DownloadError(f"Could not get smiles for {zinc_id}. Status code: {res.status_code}")
        
        info = json.loads(res.content)
        return info["smiles"]
    
    def get_catalogs(self):
        """ Get a dictionary with information about all the catalogs from ZINC.

            Returns
            -------
            dict
                Dictionary with the catalog info.
        """
        url = self._catalog_url + ".json?count=all"
        res = requests.get(url, allow_redirects=True)
        if res.status_code != requests.codes.ok:
            raise DownloadError(f"Could not download catalogs. Status code: {res.status_code}")
        
        return json.loads(res.content)
    
    def _download_zinc_file(self, file_name, url):
        """ Download a file from ZINC.
        
            Parameters
            ----------
            file_name : str
                Name that will be given to the file.
                
            url : str
                THe url from which the file will be downloaded.
        """
        res = requests.get(url, allow_redirects=True)
        if res.status_code != requests.codes.ok:
            if res.status_code == 404:
                raise ZincNotFoundError(f"Couldn't download file from {url}. Status code 404")
            elif res.status_code == 408:
                raise ZincTimeoutError(f"ZINC is taking to long to respond.")
            raise DownloadError(f"Could not download file from {url}\nStatus Code: {res.status_code}")
        
        with open(file_name, "wb") as fh:
            fh.write(res.content)
            
    def _download_batch_of_files(self, urls, download_path, fileformat, url_type, tree=True, ignore_failures=True):
        """ Download a set of files from ZINC.

            urls : list of str
                List with urls
                
            download_path : str
                The path where the files will be stored
                
            fileformat : str
                The format of the files
            
            url_type : {"2D", "3D", "CS"}
                The type of the url. Needed to extract the file name.
            
            tree : bool
                Whether to use a tree directory structure or to download all the 
                files to a single folder.
                
            ignore_failures : bool
                Whether to raise an exception if a file could not be downloaded. 
        """
        # Check input arguments
        if not isinstance(download_path, str):
            raise OpenPharmacophoreTypeError(f"Incorrect type for download_path: {type(download_path)}. Expected str")
        if not isinstance(fileformat, str):
            raise OpenPharmacophoreTypeError(f"Incorrect type for fileformat: {type(fileformat)}. Expected str")
        if not isinstance(url_type, str):
            raise OpenPharmacophoreTypeError(f"Incorrect type for url_type: {type(url_type)}. Expected str")
        if not isinstance(tree, bool):
            raise OpenPharmacophoreTypeError(f"Incorrect type for tree: {type(tree)}. Expected bool")
        if not isinstance(ignore_failures, bool):
            raise OpenPharmacophoreTypeError(f"Incorrect type for ignore_failures: {type(ignore_failures)}. Expected bool")
        
        
        file_name_extractors = {
            # For urls of this kind http://files.docking.org/2D/BA/BAAA.smi
            "2D": lambda url, fformat : url[-8:],
            # For urls of this type http://files.docking.org/3D/AA/AAML/AAAAML.xaa.db2.gz
            "3D": lambda url, fformat : url.split("/")[-1],
            # Urls of this kind https://zinc.docking.org/substances/subsets/for-sale.mol2?count=all&tranche_name=AAAA
            "CS": lambda url, fformat : url[-4:] + "." + fformat
        }
        
        try:
            file_name_extractor = file_name_extractors[url_type]
        except KeyError:
            raise InvaludUrlTypeError(f"{url_type} is not a valid url type.")
        
        for url in tqdm(urls):
            
            file_name = file_name_extractor(url, fileformat)
            if tree:
                # Create a folder for each tranch
                tranch = file_name[0:2]
                subdir = os.path.join(download_path, tranch)
                if not os.path.isdir(subdir):
                    os.mkdir(subdir)
                file_path = os.path.join(subdir, file_name)
            else: 
                file_path = os.path.join(download_path, file_name)
            
            if ignore_failures:
                try:
                    self._download_zinc_file(file_path, url)
                except (ZincTimeoutError, ZincNotFoundError):
                    continue
            else:
                self._download_zinc_file(file_path, url)
    
    def download_catalog(self, file_name, catalog_name, 
                         count=1000, availability=None, bioactive=None, 
                         biogenic=None, reactivity=None):
        """ Download a catalog from ZINC.
        
            Parameters
            ----------
            file_name : str
                The name of the file where the molecules will be stored.
            
            catalog_name : str
                Name of the catalog. Use catalogs attribute to see the available catalogs 
                
            count : int or "all"
                Number of molecules that the file will contain. Using "all" is not recommended
                as it may result in a timeout error.
            
            availability : str, default is None
                The availability of the molecules.
            
            bioactive : str, default is None
                Subset of bioactivity and drugs.
                
            biogenic : str, default is None
                Subset of biogenic.
                
            reactivity: str, default is None
                The reactivity of the molecules. 
            
        """
        url = self._get_catalog_url(file_name, catalog_name, count, availability, bioactive, biogenic, reactivity)
        self._download_zinc_file(file_name, url) 
       
        
    def download_substances(self, file_name, count=1000, availability=None, 
                         bioactive=None, biogenic=None, reactivity=None):
        """ Download a set of substances from zinc.
        
            Parameters
            ----------
            file_name : str
                The name of the file where the molecules will be stored.
                
            count : int or "all"
                Number of molecules that the file will contain. Using "all" is not recommended
                as it may result in a timeout error.
            
            availability : str, default is None
                The availability of the molecules.
            
            bioactive : str, default is None
                Subset of bioactivity and drugs.
                
            biogenic : str, default is None
                Subset of biogenic.
                
            reactivity: str, default is None
                The reactivity of the molecules. 
            
        """
        fileformat = file_name.split(".")[-1]
        url = self._append_filters_to_url(self._substances_url, fileformat, 
                                                  count, availability, bioactive, 
                                                  biogenic, reactivity)
        self._download_zinc_file(file_name, url)
        
    def download_predifined_subset(self, download_path, subset, fileformat, tree=True, ignore_failures=True):
        """ Download one of ZINC's predifined subsets.
        
            Predifined substs can only be downloaded in the following formats: smi, txt, sdf, mol2 and db2.
        
            Parameters
            ----------
            download_path : str
                The path were files will be downloaded.
            
            subset : str
                Name of the subset
            
            fileformat : str
                The format of the files that will be downloaded. Use mol2 or sdf for 3D molecules, otherwise
                use any of the other formats for 2D molecules (smiles).
            
            tree : bool
                Whether to use a tree directory structure or to download all the 
                files to a single folder.
                
            ignore_failures : bool
                Whether to raise an exception if a file could not be downloaded. 
        """
        col_list, row_list = self._predefined_subset_tranches(subset)
        
        formats_2d = ["smi", "txt"]
        formats_3d = ["sdf", "mol2", "db2"] 
        
        if fileformat in formats_2d:
            url_list = self._urls_for_tranches_2d(col_list, row_list, fileformat)
            self._download_batch_of_files(url_list, download_path, fileformat, "2D",  tree, ignore_failures)
        elif fileformat in formats_3d:
            url_list = self._urls_for_tranches_3d(col_list, row_list, fileformat)
            self._download_batch_of_files(url_list, download_path, fileformat, "3D",  tree, ignore_failures)
        else:
            raise InvalidFileFormat(f"{fileformat} is not a valid fileformat")
        
    def download_custom_subset(self, download_path, fileformat, mw_range, logp_range, tree=True, ignore_failures=True,
                                  availability=None, bioactive=None, biogenic=None, reactivity=None):
        """ Download subset with a custom molecular weight range and logP range from ZINC.
        
            This method accepts all file formats as specified in the attribute fileformats.
        
            Parameters
            ----------
            download_path : str
                The path were files will be downloaded.
                
            mw_range : 2-tuple of float
                Range of molecular weight in daltons for the downloaded molecules.
        
            logp_range : 2-tuple of float
                Range of logP for the downloaded molecules.
            
            tree : bool
                Whether to use a tree directory structure or to download all the 
                files to a single folder.
                
            ignore_failures : bool
                Whether to raise an exception if a file could not be downloaded. 
                
            availability : str, default is None
                The availability of the molecules.
            
            bioactive : str, default is None
                Subset of bioactivity and drugs.
                
            biogenic : str, default is None
                Subset of biogenic.
                
            reactivity: str, default is None
                The reactivity of the molecules.      
        """
        
        formats_2d = ["xml" ,"csv","js","json","db","solv"]
        formats_3d = ["sdf", "mol2", "db2"] 
        
        col_list, row_list = self._mw_and_logp_tranches(mw_range, logp_range)
        
        if any([availability, bioactive, biogenic, reactivity]):
            url_list = self._tranche_with_filters_url_list(col_list, row_list, availability, 
                                                        bioactive, biogenic, reactivity, 
                                                        fileformat=fileformat)
            self._download_batch_of_files(url_list, download_path, fileformat, "CS", tree, ignore_failures) 
        else:
            if fileformat == "smi" or fileformat == "txt":
                url_list = self._urls_for_tranches_2d(col_list, row_list, fileformat)
                self._download_batch_of_files(url_list, download_path, fileformat, "2D", tree, ignore_failures)    
            elif fileformat in formats_2d:
                url_list = self._tranche_with_filters_url_list(col_list, row_list, availability, 
                                                        bioactive, biogenic, reactivity, 
                                                        fileformat=fileformat)
                self._download_batch_of_files(url_list, download_path, fileformat, "CS", tree, ignore_failures)
            elif fileformat in formats_3d:
                url_list = self._urls_for_tranches_3d(col_list, row_list, fileformat)
                self._download_batch_of_files(url_list, download_path, fileformat, "3D", tree, ignore_failures)
            else:
                raise InvalidFileFormat(f"{fileformat} is not a valid file format.")
        
    ### Methods for generating urls ###
    def _append_filters_to_url(self, url, fileformat, count=1000, availability=None, 
                         bioactive=None, biogenic=None, reactivity=None):
        """ Modifies an url by adding the specified filters, count and file format.
        
            Parameters
            ----------
            url : str
                The base url
            
            fileformat : str
                The format of the file that will be downloaded
                
            count : int or "all"
                Number of molecules that the file will contain. Using "all" is not recommended
                as it may result in a timeout error.
            
            availability : str, default is None
                The availability of the molecules.
            
            bioactive : str, default is None
                Subset of bioactivity and drugs.
                
            biogenic : str, default is None
                Subset of biogenic.
                
            reactivity: str, default is None
                The reactivity of the molecules. 
            
            Returns
            -------
            url : str
                The modified url with the fileformat, count and filters.
            
        """
        if count != "all":
            if not isinstance(count, int):
                raise CountTypeError
            if count < 0:
                raise NegativeCountError
        
        self._validate_filters(fileformat, availability, bioactive, biogenic, reactivity)
        
        # Check if there are any filters
        filters = []
        if availability is not None:
            filters.append(availability)
        if bioactive is not None:
            filters.append(bioactive)
        if biogenic is not None:
            filters.append(biogenic)
        if reactivity is not None:
            filters.append(reactivity)
            
        if len(filters) > 1:
            url += "/subsets/"
            for ii in range(len(filters) - 1):
                url += filters[ii] + "+"
            url += filters[-1]       
        elif len(filters) == 1:
            url += "/subsets/" + filters[0]
        
        url += f".{fileformat}?count={count}"
        return url
    
    def _get_catalog_url(self, file_name, catalog_name, count=1000, 
                         availability=None, bioactive=None, 
                         biogenic=None, reactivity=None):
        """ Obtain the url for downloading a catalog."""
        
        if catalog_name not in self.catalogs.keys():
            raise InvalidCatalogError(f"{catalog_name} is not a valid catalog name")
        
        catalog_short_name = self.catalogs[catalog_name]
        
        fileformat = file_name.split(".")[-1]
        
        url = self._catalog_url + "/" + catalog_short_name + "/substances"
        url = self._append_filters_to_url(url, fileformat, count, availability, 
                                                bioactive, biogenic, reactivity)
        
        return url
    
    def _urls_for_tranches_2d(self, col_list, row_list, fileformat="smi"):
        """ Returns a list of urls to download files in smi format for the specified tranches.
        
            Parameters
            ----------
            col_list : list of str
                List with the columns names. Columns are named with letters from A
                to K. They correspond to molecular weight.
                
            row_list : list of str
                List with the row names. Rows are named with letters from A to K.
                They correspond to LogP.
                
            Returns
            -------
            url_list : list of str
                The urls. 
            
        """
        url_list = []
        if fileformat != "smi" and fileformat != "txt":
            raise InvalidFileFormat(f"{fileformat} is not a valid file format. Valid formats are smi or txt")
        
        for col in col_list:
            for row in row_list:           
                tranch = col + row
                # Each tranch is divided into various files from A to E
                tranch_subcategories = itertools.product("ABCE", "ABCD", repeat=1)
                url = self._tranches_2d_url + tranch + "/" + tranch      
                for subtranch in tranch_subcategories:
                    url_download = url + subtranch[0] + subtranch[1] + "." + fileformat
                    url_list.append(url_download)

        return url_list
    
    def _urls_for_tranches_3d(self, col_list, row_list, fileformat):
        """ Get a list of urls to download files in a 3D format for the specified tranches.
        
            Parameters
            ----------
            col_list : list of str
                List with the columns names. Columns are named with letters from A
                to K. They correspond to molecular weight.
                
            row_list : list of str
                List with the row names. Rows are named with letters from A to K.
                They correspond to LogP.
        
            fileformat : {"sdf", "mol2", "db2"}
                The format of the files.
                
            Returns
            -------
            url_list : list of str
                The urls. 
        
        """
        formats_3d = ["sdf", "mol2", "db2"] 
        if fileformat not in formats_3d:
            raise InvalidFileFormat(f"{fileformat} is not a valid 3D file format. Valid formats are: {formats_3d}")
        
        tranches = []
        for column in col_list:
            for row in row_list:
                tranches.append(column + row)
                
        base_url = "http://files.docking.org/3D/"
        
        urls3d_dir = "./data/zinc/urls3d/"
        url_list = []
        for tranch in tranches:
            file = pkg_resources.resource_filename("openpharmacophore", urls3d_dir + tranch + ".uri")
            with open(file, "r") as fh:
                for line in fh.readlines():
                    url = base_url + tranch + "/" + line.rstrip() + "." + fileformat + ".gz"
                    url_list.append(url)
        
        return url_list
    
    def _tranch_url_with_filters(self, tranch, availability, bioactive, biogenic, reactivity, fileformat="smi"):
        """ Generate a url for a tranch with filters.
        """
        
        url = self._append_filters_to_url(self._substances_url, fileformat=fileformat,
                                            count="all", availability=availability,
                                            bioactive=bioactive, biogenic=biogenic,
                                            reactivity=reactivity)
        url += f"&tranche_name={tranch}"        
      
        return url
    
    def _tranche_with_filters_url_list(self, col_list, row_list, availability, bioactive, biogenic, reactivity, fileformat):
        """ Generate a list of urls for a set of tranches with filters.
        """
        # Get urls for smi files
        url_list = []
        for col in col_list:
            for row in row_list:
                tranch = col + row     
                tranch_subcategories = itertools.product("ABCE", "ABCD", repeat=1)
                for subtranch in tranch_subcategories:
                    url = self._tranch_url_with_filters(tranch + subtranch[0] + subtranch[1], 
                                                        availability, 
                                                        bioactive, 
                                                        biogenic, 
                                                        reactivity,
                                                        fileformat)
                    url_list.append(url)
                        
        return url_list
    
    ### Methods for generating tranches ###

    def _predefined_subset_tranches(self, subset):
        """ Obtain the tranches for one of ZINC's predifined subsets.

            Parameters
            ----------
            subset : str
                Name of the subset.
            
            Returns
            -------
            col_list : list of str
                List with the tranches columns names. Columns are named with letters from A
                to K. They correspond to molecular weight.
                
            row_list : list of str
                List with the tranches row names. Rows are named with letters from A to K.
                They correspond to LogP.
            
        """
        # Subsets defined by ZINC
        subsets = {
            # First tuple is start and end columns, second tuple is start and end rows
            "Drug-Like"   : [(1, 9), (0, 9)],
            "Lead-Like"   : [(2, 4), (0, 6)],
            "Lugs"        : [(4, 8), (0, 8)],
            "Goldilocks"  : [(2, 4), (3, 5)],
            "Fragments"   : [(0, 1), (0, 6)],
            "Flagments"   : [(1, 3), (0, 6)],
            "Big-n-Greasy": [(9, 10), (8, 10)],
            "Shards"      : [(0, 0), [0, 10]],
        }
        
        if subset not in list(subsets.keys()):
            raise InvalidSubsetError(f"{subset} is not a valid subset. Valid subsets are: {list(subsets.keys())}")
        
        # ZINC molecular weight is categorized in columns from A to K
        mw_cols = list(string.ascii_uppercase[:11])
        # LogP is categorized in rows from A to K
        logp_rows = mw_cols
        
        start_col, end_col = subsets[subset][0]
        start_row, end_row = subsets[subset][1]
        
        col_list = mw_cols[start_col:end_col + 1]
        row_list = logp_rows[start_row:end_row + 1]

        return col_list, row_list
                        
    def _mw_and_logp_tranches(self, mw_range, logp_range):
        """ Obtain the tranches for a custom subset with the specified molecular weight and logP
            range.

            Parameters
            ----------
            mw_range : 2-tuple of float
            Range of molecular weight in daltons for the downloaded molecules.
        
            logp_range : 2-tuple of float
                Range of logP for the downloaded molecules.
            
            Returns
            -------
            col_list : list of str
                List with the tranches columns names. Columns are named with letters from A
                to K. They correspond to molecular weight.
                
            row_list : list of str
                List with the tranches row names. Rows are named with letters from A to K.
                They correspond to LogP.
            
        """
        # This are the values that ZINC accepts
        mw_values = [200, 250, 300, 325, 350, 375, 400, 425, 450, 500, 550]
        logp_values = [-1, 0, 1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6]
        
        # Validate molecular weight
        if mw_range[0] > mw_range[1]:
            raise InvalidMolecularWeightRangeError("First number must be smaller")
        if mw_range[0] < mw_values[0] or mw_range[0] > mw_values[-1]:
            raise  InvalidMolecularWeightRangeError(f"Molecular weight must be a value between {mw_values[0]} and {mw_values[-1]}")
        if mw_range[1] < mw_values[0] or mw_range[1] > mw_values[-1]:
            raise  InvalidMolecularWeightRangeError(f"Molecular weight must be a value between {mw_values[0]} and {mw_values[-1]}")
        
        # Validate logP
        if logp_range[0] > logp_range[1]:
            raise InvalidLogPRangeError("First number must be smaller")
        if logp_range[0] < logp_values[0] or logp_range[0] > logp_values[-1]:
            raise  InvalidLogPRangeError(f"LogP must be a value between {logp_values[0]} and {logp_values[-1]}")
        if logp_range[1] < logp_values[0] or logp_range[1] > logp_values[-1]:
            raise  InvalidLogPRangeError(f"LogP must be a value between {logp_values[0]} and {logp_values[-1]}")
        
        # ZINC molecular weight is categorized in columns from A to K
        mw_cols = list(string.ascii_uppercase[:11])
        # LogP is categorized in rows from A to K
        logp_rows = mw_cols
        
        mw_lower_bound, mw_upper_bound = mw_range
        logp_lower_bound, logp_upper_bound = logp_range
        
        molecular_weight = dict(zip(mw_values, mw_cols))
        logp = dict(zip(logp_values, logp_rows))
        
        # Discretize mol weight and logp values:
        mw_lower_bound = self.discretize_values(mw_lower_bound, mw_values, "Molecular weight")
        mw_upper_bound = self.discretize_values(mw_upper_bound, mw_values, "Molecular weight", lower=False)
        logp_lower_bound = self.discretize_values(logp_lower_bound, logp_values, "LogP")
        logp_upper_bound = self.discretize_values(logp_upper_bound, logp_values, "LogP", lower=False)

        start_col = mw_cols.index(molecular_weight[mw_lower_bound])
        end_col = mw_cols.index(molecular_weight[mw_upper_bound]) 
    
        start_row = logp_rows.index(logp[logp_lower_bound])
        end_row = logp_rows.index(logp[logp_upper_bound])
        
        col_list = mw_cols[start_col:end_col + 1]
        row_list = logp_rows[start_row:end_row + 1]
        
        return col_list, row_list
    
    ### Methods to obtain data necessary for the class to work (i.e catalog names, filter names)
    def _get_filters(self):
        """Get all the filters available in ZINC.
        
           Returns
           -------
           dict
            Dictionary of filters grouped by category
        """
        filters = {
            "Availability"      : ["bb","for-sale","not-for-sale","now","wait-ok"],
            "BioactiveAndDrugs" : ["fda","in-cells","in-trials","in-vivo","in-vitro","world"],
            "Biogenic"          : ["biogenic","endogenous","metabolites"],
            "Reactivity"        : ["anodyne","clean","hot-ok","reactive-ok","standard-ok"]
        }
        
        return filters

    def _get_catalogs_names(self):
        """ Return the names and short names of all catalogs. 
            Short names are used in the url for searching.
            
            Returns
            -------
            dict
                Dictionary with catalog full name as keys and short name as values.
        """
        return {cat["name"] : cat["short_name"] for cat in self.get_catalogs()}
    
    def _load_catalogs_from_file(self):
        """ Load the catalogs from a file in case they can't be downloaded"""
        catalogs_file = pkg_resources.resource_filename("openpharmacophore", "./data/zinc/catalogs.json")
        with open(catalogs_file, "r") as fh:
            catalogs = json.load(fh)
        return catalogs
        
    ### Auxiliary Methods ###

    @staticmethod
    def discretize_values(value, bins, name, lower=True):
        """Discretize a molecualr weight or logp value.

        Parameters
        ----------
        value : int
            Value that will be discretized.

        bins : list of int
            List containing the bins to discretize the value
        
        name : str
            Name of the variable that will be discretized
        
        lower : bool
            If True the lower bound will be assigned to value.
            Else the upper bound will be assigned.

        Returns
        -------
        value : double
            The discretized value.

        """
        for ii in range(len(bins) - 1):
            if value < bins[0]:
                raise OpenPharmacophoreValueError("{} must be at least {}".format(name, bins[0]))
            elif value >= bins[-1]:
                value = bins[-1]
            if value > bins[ii] and value < bins[ii + 1]:
                if lower:
                    return bins[ii]
                else:
                    return bins[ii + 1]

        return value
    
    def _validate_filters(self, fileformat, availability=None, bioactive=None, biogenic=None, reactivity=None):
        """ Validate the filters passed to the download_substances and download_catalogs methods.
        
            Parameters
            ----------
            availability : str, default is None
                The availability of the molecules.
            
            bioactive : str, default is None
                Subset of bioactivity and drugs.
                
            biogenic : str, default is None
                Subset of biogenic.
                
            reactivity: str, default is None
                The reactivity of the molecules. 
                
            Raises
            ------
            InvalidFileFormat
            
            InvalidAvailabilityError
            
            InvalidBioactiveError
            
            InvalidBiogenicError
            
            InvalidReactivityError
        """
        if fileformat not in self.file_formats:
            raise InvalidFileFormat(f"{fileformat} is not a valid fileformat.")
        
        if availability:
            if availability not in self.filters["Availability"]:
                raise InvalidAvailabilityError(f"{availability} is not a valid availability.")
        
        if bioactive:
            if bioactive not in self.filters["BioactiveAndDrugs"]:
                raise InvalidBioactiveError(f"{bioactive} is not a valid bioactivity.")
        
        if biogenic:    
            if biogenic not in self.filters["Biogenic"]:
                raise InvalidBiogenicError(f"{biogenic} is not a valid biogenic.")
        
        if reactivity:    
            if reactivity not in self.filters["Reactivity"]:
                raise InvalidReactivityError(f"{reactivity} is not a valid reactivity.")
