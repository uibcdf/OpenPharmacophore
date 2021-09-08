import os
import requests
import string

def download_ZINC2D_smiles(download_path, subset="Drug-like", mol_weight_range=None, logp_range=None, **kwargs):
    """Download a set of molecules as .smi files from ZINC database.

        Parameters
        ----------
        download_path: str
            Name of the path where the files will be downloaded.

        subset: str (optional)
            Name of the ZINC subset to be downloaded. Available subsets are "Drug-Like",
            "Lead-Like", "Lugs", "Goldilocks", "Fragments", "Flagments", "Big-n-Greasy"
            and "Shards".
        
        mol_weight_range: 2-tuple of int
            Range of molecular weight for the downloaded molecules.
        
        logp_range: 2-tuple of int
            Range of logP for the downloaded molecules.

        Note
        -------
        Nothing is returned. Files are downloaded

        """
    if (mol_weight_range is None or logp_range is None) and subset is None:
        raise ValueError("Missing parameters")
    
    # Subsets defined by ZINC
    subsets = {
        # First tuple is start and end columns, second tuple is start and end rows
        "Drug-Like": [(1, 9), (0, 9)],
        "Lead-Like": [(2, 4), (0, 7)],
        "Lugs": [(4, 8), (0, 7)],
        "Goldilocks": [(2, 4), (3, 5)],
        "Fragments": [(0, 1), (0, 6)],
        "Flagments": [(1, 3), (0, 6)],
        "Big-n-Greasy": [(9, 10), (8, 10)],
        "Shards": [(0, 0), [0, 10]]
    }

    # This are the values that ZINC accepts
    mw_values = [200, 250, 300, 325, 350, 375, 400, 425, 450, 500, 550]
    logp_values = [-1, 0, 1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6]
    # ZINC molecular weight is categorized in columns from A to K
    mw_cols = list(string.ascii_uppercase[:11])
    # LogP is categorized in rows from A to K
    logp_rows = mw_cols
    
    molecular_weight = dict(zip(mw_values, mw_cols))
    logp = dict(zip(logp_values, logp_rows))
    
    if subset is None:
        mw_lower_bound, mw_upper_bound = mol_weight_range
        logp_lower_bound, logp_upper_bound = logp_range
        
        # Discretize mol weight and logp values:
        mw_lower_bound = discretize_values(mw_lower_bound, mw_values, "Molecular weight")
        mw_upper_bound = discretize_values(mw_upper_bound, mw_values, "Molecular weight", lower=False)
        logp_lower_bound = discretize_values(logp_lower_bound, logp_values, "LogP")
        logp_upper_bound = discretize_values(logp_upper_bound, logp_values, "LogP", lower=False)

        start_col = mw_cols.index(molecular_weight[mw_lower_bound])
        end_col = mw_cols.index(molecular_weight[mw_upper_bound]) 
       
        start_row = logp_rows.index(logp[logp_lower_bound])
        end_row = logp_rows.index(logp[logp_upper_bound])
    else:
        start_col, end_col = subsets[subset][0]
        start_row, end_row = subsets[subset][1]
    
    col_list = mw_cols[start_col:end_col + 1]
    row_list = logp_rows[start_row:end_row + 1]

    # For testing purposes
    if kwargs:
        if kwargs["testing"] == True:
            url_list = []
    
    base_url = "http://files.docking.org/2D/"
    # Get urls and download files
    for col in col_list:
        for row in row_list:
            tranch = col + row
            url = base_url + tranch + "/" + tranch
            # Each tranch is divided into various files from A to E
            for f in ["A", "B", "C", "E"]:
                for j in ["A", "B"]:
                    url_download = url + f + j + ".smi"

                    if kwargs:
                        if kwargs["testing"] == True:
                            url_list.append(url_download)
                            continue
                    try:
                        r = requests.get(url_download, allow_redirects=True)
                    except:
                        print("Could not download file from {}".format(url_download))
                    file_name = tranch + f + j + ".smi"
                    file_path = os.path.join(download_path, file_name)
                    with open(file_path, "wb") as file:
                        file.write(r.content)

    if kwargs:
        if kwargs["testing"] == True:
            return url_list


def discretize_values(value, bins, name, lower=True):
    """Download a set of molecules as .smi files from ZINC database.

    Parameters
    ----------
    value: int
        Value that will be disctretized.

    bins: list of int
        List containing the bins ti discretize the value
    
    name: str
        Name of the variable that will be discretized
    
    lower: bool
        If True the lower bound will be assigned to value.
        Else the upper bound will be assigned.

    Returns
    -------
    Nothing is returned. Files are downloaded

    """
    for i in range(len(bins) - 1):
        if value < bins[0]:
            raise ValueError("{} must be at least {}".format(name, bins[0]))
        elif value >= bins[-1]:
            value = bins[-1]
        if value > bins[i] and value < bins[i + 1]:
            if lower:
                value = bins[i]
            else:
                value = bins[i + 1]

    return value