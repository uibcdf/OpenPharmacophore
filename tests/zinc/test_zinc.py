import openpharmacophore.zinc.zinc as zinc
import openpharmacophore._private_tools.exceptions as exc
import pytest
import os


def test_download_file_success(mocker):
    mock_get = mocker.patch("openpharmacophore.zinc.zinc.requests.get")
    mock_open = mocker.patch("openpharmacophore.zinc.zinc.open", new=mocker.mock_open())

    mock_get.return_value.status_code = 200
    mock_get.return_value.content = "c1ccccc1"

    zinc._download_file("test_file.smi", "http://files.docking.org/3D/BA/AAML/BAAAML.smi")
    mock_open.assert_called_with("test_file.smi", "wb")
    mock_open.return_value.write.assert_called_once_with("c1ccccc1")


def test_download_file_unsuccessful_raises_error(mocker):
    mock_get = mocker.patch("openpharmacophore.zinc.zinc.requests.get")
    mock_get.return_value.status_code = 404
    with pytest.raises(exc.ZincDownloadError):
        zinc._download_file("test_file.smi", "http://files.docking.org/3D/BA/AAML/BAAAML.smi")


def test_download_multiple_files_3d_url(mocker):
    mock_download = mocker.patch("openpharmacophore.zinc.zinc._download_file")
    mock_mkdir = mocker.patch("openpharmacophore.zinc.zinc.os.mkdir")

    urls = [
        "http://files.docking.org/3D/BA/AAML/BAAAML.sdf",
        "http://files.docking.org/3D/BB/AAML/BBAAMM.sdf",
        "http://files.docking.org/3D/CA/AAML/CAAAMN.sdf"
    ]
    download_path = "./data"
    zinc._download_multiple_files(urls, "./data")

    arg_list = mock_mkdir.call_args_list
    assert len(arg_list) == 3
    assert arg_list[0] == mocker.call(os.path.join(download_path, "BA"))
    assert arg_list[1] == mocker.call(os.path.join(download_path, "BB"))
    assert arg_list[2] == mocker.call(os.path.join(download_path, "CA"))

    arg_list = mock_download.call_args_list
    assert len(arg_list) == 3
    assert arg_list[0] == mocker.call(os.path.join(download_path, "BA", "BAAAML.sdf"), urls[0])
    assert arg_list[1] == mocker.call(os.path.join(download_path, "BB", "BBAAMM.sdf"), urls[1])
    assert arg_list[2] == mocker.call(os.path.join(download_path, "CA", "CAAAMN.sdf"), urls[2])


def test_download_multiple_files_2d_url(mocker):
    mock_download = mocker.patch("openpharmacophore.zinc.zinc._download_file")
    mock_mkdir = mocker.patch("openpharmacophore.zinc.zinc.os.mkdir")

    urls = [
        "http://files.docking.org/2D/DD/DDAA.smi",
        "http://files.docking.org/2D/ED/EDAA.smi",
        "http://files.docking.org/2D/FD/FDAA.smi"
    ]
    download_path = "./data"
    zinc._download_multiple_files(urls, "./data")

    arg_list = mock_mkdir.call_args_list
    assert len(arg_list) == 3
    assert arg_list[0] == mocker.call(os.path.join(download_path, "DD"))
    assert arg_list[1] == mocker.call(os.path.join(download_path, "ED"))
    assert arg_list[2] == mocker.call(os.path.join(download_path, "FD"))

    arg_list = mock_download.call_args_list
    assert len(arg_list) == 3
    assert arg_list[0] == mocker.call(os.path.join(download_path, "DD", "DDAA.smi"), urls[0])
    assert arg_list[1] == mocker.call(os.path.join(download_path, "ED", "EDAA.smi"), urls[1])
    assert arg_list[2] == mocker.call(os.path.join(download_path, "FD", "FDAA.smi"), urls[2])


@pytest.fixture()
def mock_tranche_dict_2d():
    return {
        "AA": ["AAAA", "AAAB"],
        "BA": ["BAAA", "BAAB"],
        "CA": ["CAAA", "CAAB"],
        "DA": ["DAAA", "DAAB"],
        "AB": ["ABAA", "ABAB"],
        "BB": ["BBAA", "BBAB"],
        "CB": ["CBAA", "CBAB"],
        "DB": ["DBAA", "DBAB"],
    }


@pytest.fixture()
def mock_tranche_dict_3d():
    return {
        "AA": ["AAAAML.xaa", "AAAAMM.xaa"],
        "BA": ["BAAAML.xaa", "BAAAMM.xaa"],
        "CA": ["CAAAML.xaa", "CAAAMM.xaa"],
        "DA": ["DAAAML.xaa", "DAAAMM.xaa"],
        "AB": ["ABAAML.xaa", "ABAAMM.xaa"],
        "BB": ["BBAAML.xaa", "BBAAMM.xaa"],
        "CB": ["CBAAML.xaa", "CBAAMM.xaa"],
        "DB": ["DBAAML.xaa", "DBAAMM.xaa"],
    }


def test_load_tranches_dict_2d():
    tranches_2d = zinc._load_tranches_dict("2d")
    assert len(tranches_2d) == 121
    assert "AA" in tranches_2d
    assert "JC" in tranches_2d
    assert len(tranches_2d["AK"]) == 4


def test_load_tranches_dict_3d():
    tranches_3d = zinc._load_tranches_dict("3d")
    assert len(tranches_3d) == 121
    assert "AA" in tranches_3d
    assert "JC" in tranches_3d
    assert len(tranches_3d["AJ"]) == 4


def test_closest_value():
    values = [200, 250, 300, 325, 350, 375, 400, 425, 450, 500, 550]
    assert zinc._closest_value(225, values) == 200
    assert zinc._closest_value(360, values) == 350
    assert zinc._closest_value(100, values) == 200
    assert zinc._closest_value(600, values) == 550


def test_tranche_rows_and_cols():
    rows, cols = zinc._tranches_rows_and_cols(
        logp=(-1, 0), mw=(250, 300))
    assert rows == (0, 1)
    assert cols == (1, 2)


def test_urls_from_mw_and_logp_2d(mock_tranche_dict_2d):
    urls = zinc._urls_from_mw_and_logp(
        logp=(-1, 0), mw=(250, 300), tranches=mock_tranche_dict_2d)
    assert urls == [
        "BA/BAAA",
        "BA/BAAB",
        "CA/CAAA",
        "CA/CAAB",
        "BB/BBAA",
        "BB/BBAB",
        "CB/CBAA",
        "CB/CBAB",
    ]


def test_urls_from_mw_and_logp_3d(mock_tranche_dict_3d):
    urls = zinc._urls_from_mw_and_logp(
        logp=(-1, 0), mw=(250, 300), tranches=mock_tranche_dict_3d
    )
    assert urls == [
        "BA/BAAAML.xaa",
        "BA/BAAAMM.xaa",
        "CA/CAAAML.xaa",
        "CA/CAAAMM.xaa",
        "BB/BBAAML.xaa",
        "BB/BBAAMM.xaa",
        "CB/CBAAML.xaa",
        "CB/CBAAMM.xaa",
    ]
