from openpharmacophore import StructureBasedPharmacophore, PharmacophoricPoint
import openpharmacophore.data as data
import nglview as nv
import numpy as np
import pyunitwizard as puw
import pytest
from copy import deepcopy
import os


def test_init_from_pdb_file():
    pharmacophore = StructureBasedPharmacophore(data.pdb["1ncr"])
    assert pharmacophore.num_frames == 1
    assert len(pharmacophore) == 0
    assert pharmacophore._pdb == data.pdb["1ncr"]


def test_init_from_pdb_id(mocker):
    mocker.patch("openpharmacophore.StructureBasedPharmacophore._fetch_pdb",
                 return_value="pdb")
    pharmacophore = StructureBasedPharmacophore("1NCR")
    assert len(pharmacophore) == 0
    assert pharmacophore.num_frames == 1
    assert pharmacophore._pdb == "pdb"


def test_is_pdb_id():
    assert StructureBasedPharmacophore._is_pdb_id("1A52")
    assert StructureBasedPharmacophore._is_pdb_id("5UFW")
    assert not StructureBasedPharmacophore._is_pdb_id("1A52D")
    assert not StructureBasedPharmacophore._is_pdb_id("1A3")


def test_fetch_pdb(mocker):
    mock_get = mocker.patch("requests.get")
    mock_get.return_value.status_code = 200
    mock_get.return_value.content = b"pdb"

    assert StructureBasedPharmacophore._fetch_pdb("1A52") == "pdb"
    mock_get.assert_called_once_with('http://files.rcsb.org/download/1A52.pdb', allow_redirects=True)


def test_from_file_moe():
    pharmacophore = StructureBasedPharmacophore(None)
    pharmacophore.from_file(data.pharmacophores["gmp"])
    assert len(pharmacophore) == 1
    assert len(pharmacophore[0]) == 10


def test_from_file_mol2():
    pharmacophore = StructureBasedPharmacophore(None)
    pharmacophore.from_file(data.pharmacophores["elastase"])
    assert len(pharmacophore) == 8
    assert pharmacophore.num_frames == 8
    assert pharmacophore._pharmacophores_frames == list(range(0, 8))


def test_from_file_json():
    pharmacophore = StructureBasedPharmacophore(None)
    pharmacophore.from_file(data.pharmacophores["1M70"])
    assert len(pharmacophore) == 1
    assert len(pharmacophore[0]) == 5


def test_from_file_pml():
    pharmacophore = StructureBasedPharmacophore(None)
    pharmacophore.from_file(data.pharmacophores["ligscout"])
    assert len(pharmacophore) == 1
    assert len(pharmacophore[0]) == 4


@pytest.fixture
def pharmacophore_one_frame():
    acceptor = PharmacophoricPoint(
        "hb acceptor",
        puw.quantity([1.0, 1.0, 1.0], "angstroms"),
        puw.quantity(1.0, "angstroms")

    )
    donor = PharmacophoricPoint(
        "hb donor",
        puw.quantity([2.0, 0.0, 3.0], "angstroms"),
        puw.quantity(1.5, "angstroms")
    )
    aromatic = PharmacophoricPoint(
        "aromatic ring",
        puw.quantity([0.0, 1.5, 2.0], "angstroms"),
        puw.quantity(1.0, "angstroms")
    )

    ph = StructureBasedPharmacophore(None)
    ph.add_frame()
    ph.add_points_to_frame([acceptor, donor, aromatic], 0)
    return ph


def test_add_points_to_frame(pharmacophore_one_frame):
    ph = pharmacophore_one_frame

    assert ph.num_frames == 1
    assert len(ph) == 1
    assert len(ph[0]) == 3
    assert ph[0][0].feature_name == "hb acceptor"
    assert ph[0][1].feature_name == "hb donor"
    assert ph[0][2].feature_name == "aromatic ring"


def test_remove_point(pharmacophore_one_frame):
    ph = deepcopy(pharmacophore_one_frame)
    ph.remove_point(1, 0)
    assert len(ph[0]) == 2
    assert ph[0][0].feature_name == "hb acceptor"
    assert ph[0][1].feature_name == "aromatic ring"


def test_adding_single_frame_to_view_updates_components(pharmacophore_one_frame):
    view = nv.NGLWidget()
    assert len(view._ngl_component_ids) == 0
    pharmacophore_one_frame.add_to_view(view, 0)
    assert len(view._ngl_component_ids) == 3


def assert_file_is_created(file_name):
    assert os.path.isfile(file_name)
    os.remove(file_name)


def test_to_json(pharmacophore_one_frame):
    file_name = "ph.json"
    pharmacophore_one_frame.to_json(file_name, 0)
    assert_file_is_created(file_name)


def test_to_ligand_scout(pharmacophore_one_frame):
    file_name = "ph.pml"
    pharmacophore_one_frame.to_ligand_scout(file_name, 0)
    assert_file_is_created(file_name)


def test_to_moe(pharmacophore_one_frame):
    file_name = "ph.ph4"
    pharmacophore_one_frame.to_moe(file_name, 0)
    assert_file_is_created(file_name)


def test_to_mol2(pharmacophore_one_frame):
    file_name = "ph.mol2"
    pharmacophore_one_frame.to_mol2(file_name, 0)
    assert_file_is_created(file_name)


def test_to_rdkit(pharmacophore_one_frame):
    rdkit_ph, radii = pharmacophore_one_frame.to_rdkit(0)

    assert len(radii) == 3
    assert np.all(np.array(radii) == np.array([1.0, 1.5, 1.0]))

    feats = rdkit_ph.getFeatures()
    assert len(feats) == 3

    acceptor = feats[0]
    assert acceptor.GetFamily() == "Acceptor"
    assert np.allclose(acceptor.GetPos().x, 1.0)
    assert np.allclose(acceptor.GetPos().y, 1.0)
    assert np.allclose(acceptor.GetPos().z, 1.0)

    ring_1 = feats[1]
    assert ring_1.GetFamily() == "Donor"
    assert np.allclose(ring_1.GetPos().x, 2.0)
    assert np.allclose(ring_1.GetPos().y, 0.0)
    assert np.allclose(ring_1.GetPos().z, 3.0)

    ring_2 = feats[2]
    assert ring_2.GetFamily() == "Aromatic"
    assert np.allclose(ring_2.GetPos().x, 0.0)
    assert np.allclose(ring_2.GetPos().y, 1.5)
    assert np.allclose(ring_2.GetPos().z, 2.0)


@pytest.fixture()
def analyzed_pharmacophore():
    pharmacophore = StructureBasedPharmacophore(data.pdb["1ncr"])
    pharmacophore.analyze()
    return pharmacophore


def test_protein_ligand_interactions(analyzed_pharmacophore):
    interactions = analyzed_pharmacophore._interactions
    assert len(interactions) == 2
    assert "MYR:D:4000" in interactions
    assert "W11:A:7001" in interactions


def test_points_from_interactions(analyzed_pharmacophore):
    interactions = analyzed_pharmacophore._interactions
    points = StructureBasedPharmacophore._points_from_interactions(
        interactions, "W11:A:7001", 1.0)

    assert len(points) == 9
    n_hydrophobics = len([p for p in points if p.short_name == "H"])
    n_rings = len([p for p in points if p.short_name == "R"])
    assert n_hydrophobics == 8
    assert n_rings == 1


def test_extract_single_frame_without_analyzing():
    pharmacophore = StructureBasedPharmacophore(data.pdb["1ncr"])
    pharmacophore.extract("W11:A:7001", 1.0)
    assert len(pharmacophore) == 1
    assert len(pharmacophore[0]) == 9
    assert pharmacophore._pharmacophores_frames == [0]


@pytest.fixture()
def analyzed_and_extracted_pharmacophore():
    pharmacophore = StructureBasedPharmacophore(data.pdb["1ncr"])
    pharmacophore.analyze()
    pharmacophore.extract("W11:A:7001", 1.0)

    return pharmacophore


def test_extract_single_frame_after_analyzing(analyzed_and_extracted_pharmacophore):
    assert len(analyzed_and_extracted_pharmacophore) == 1
    assert len(analyzed_and_extracted_pharmacophore[0]) == 9


def test_interactions_and_pdb_are_deleted_after_extraction(analyzed_and_extracted_pharmacophore):
    with pytest.raises(AttributeError):
        analyzed_and_extracted_pharmacophore._interactions
    with pytest.raises(AttributeError):
        analyzed_and_extracted_pharmacophore._pdb
    with pytest.raises(AttributeError):
        analyzed_and_extracted_pharmacophore._pdb_is_str


def test_ligand_ids(analyzed_pharmacophore):
    assert analyzed_pharmacophore.ligand_ids() == [
        "W11:A:7001",
        "MYR:D:4000",
    ]
