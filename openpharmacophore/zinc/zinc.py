from .._private_tools.exceptions import ZincDownloadError
from ..data import zinc
import requests
import os
import string


def _download_file(file_name, url):
    """ Download a file from ZINC.

        Parameters
        ----------
        file_name : str
                Name that will be given to the file.

        url : str
            The url from which the file will be downloaded.
    """
    res = requests.get(url)
    if res.status_code != 200:
        raise ZincDownloadError(res.status_code, url)

    with open(file_name, "wb") as fp:
        fp.write(res.content)


def _download_multiple_files(urls, download_path="./"):
    """ Download multiple files from ZINC.

        Parameters
        ----------
        urls : List[str]
            List with urls of the files that will be downloaded.

        download_path : str, default="./"
            The path were the files will be downloaded.
    """
    for url in urls:
        file_name = url.split("/")[-1]
        tranch = file_name[0:2]
        subdir = os.path.join(download_path, tranch)
        if not os.path.isdir(subdir):
            os.mkdir(subdir)
        file_path = os.path.join(subdir, file_name)
        try:
            _download_file(file_path, url)
        except ZincDownloadError:
            continue


def _load_tranches_dict(name):
    """ Returns a dictionary with the names of the tranches
        for 2D urls.

        Returns
        -------
        tranches : dict[str, list[str]]
    """
    if name == "2d":
        path = zinc["tranches_2D"]
    elif name == "3d":
        path = zinc["tranches_3D"]
    else:
        raise ValueError
    tranches = {}
    with open(path) as fp:
        for line in fp.readlines():
            sub_tranches = line.split()
            tranche_name = line[0:2]
            tranches[tranche_name] = sub_tranches
    return tranches


def _closest_value(value, value_list):
    """ Returns the number that is closest to the given value
        in value list.

        Parameters
        ----------
        value : float
            The number that its closest value will be returned.

        value_list : list[float]
            An ordered list with the values.

        Returns
        -------
        float
            The closest number to 'value'.
    """
    if value < value_list[0]:
        return value_list[0]
    elif value > value_list[-1]:
        return value_list[-1]

    low = 0
    high = len(value_list) - 1
    while low <= high:
        mid = low + (high - low) // 2

        if value_list[mid] == value:
            return value
        # If value is between mid and mid + 1, we return mid
        elif value_list[mid] < value < value_list[mid + 1]:
            return value_list[mid]
        elif value_list[mid] < value:
            low = mid + 1
        else:
            high = mid - 1
    return -1


def _tranches_rows_and_cols(logp, mw):
    """ Get the index of the rows and columns corresponding to the
        given logp and mw values.

        Parameters
        ----------
        logp : tuple[float, float]
            The lowest value of logP and the greatest value of logP
            for the molecules.

        mw : tuple[float, float]
            The lowest value and greatest value of molecular weight
            for the molecules.

        Returns
        -------
        rows : tuple[int, int]
            The index of the first and final row.

        columns : tuple[int, int]
            The index of the first and final column.

    """
    mw_values = [200, 250, 300, 325, 350, 375, 400, 425, 450, 500, 550]
    logp_values = [-1, 0, 1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6]

    start_mw = _closest_value(mw[0], mw_values)
    end_mw = _closest_value(mw[1], mw_values)
    start_logp = _closest_value(logp[0], logp_values)
    end_logp = _closest_value(logp[1], logp_values)

    columns = (mw_values.index(start_mw), mw_values.index(end_mw))
    rows = (logp_values.index(start_logp), logp_values.index(end_logp))

    return rows, columns


def _urls_from_mw_and_logp(logp, mw, tranches):
    """ Returns a list of urls for the files containing molecules
        between the given logP and molecular weight.

        Parameters
        ----------
        logp : tuple[float, float]
            The lowest value of logP and the greatest value of logP
            for the molecules.

        mw : tuple[float, float]
            The lowest value and greatest value of molecular weight
            for the molecules.

        tranches : dict[str, list[str]]
            Dictionary with a list of the names of each sub-tranche.

        Returns
        -------
        urls : list[str]
            List with urls
    """
    urls = []

    rows, columns = _tranches_rows_and_cols(logp, mw)
    # Rows and columns are named from A to K
    row_names = string.ascii_uppercase[:11]

    for ii in range(rows[0], rows[1] + 1):
        for jj in range(columns[0], columns[1] + 1):
            tranche_name = row_names[jj] + row_names[ii]
            sub_tranches = tranches[tranche_name]
            for sub_tranche in sub_tranches:
                url = tranche_name + "/" + sub_tranche
                urls.append(url)

    return urls


def _predefined_subset_urls(subset):
    """ Get a list of urls for one of ZINC's predifined subsets.

        Parameters
        ----------
        subset : str
            Name of the subset

        Returns
        -------
        list[str]
            List with urls
    """
    subsets_dict = {
        # First tuple is logp and second column is mw
        "Drug-Like": [(-1, 5), (250, 500)],
        "Lead-Like": [(-1, 3.5), (300, 350)],
        "Lugs": [(-1, 4.5), (350, 450)],
        "Goldilocks": [(2, 3), (300, 350)],
        "Fragments": [(-1, 3.5), (200, 250)],
        "Flagments": [(-1, 3.5), (250, 325)],
        "Big-n-Greasy": [(4.5, 6), (500, 550)],
        "Shards": [(-1, 6), (200, 200)],
    }
    logp, mw = subsets_dict[subset]
    return _urls_from_mw_and_logp(logp, mw, _load_tranches_dict("2d"))
