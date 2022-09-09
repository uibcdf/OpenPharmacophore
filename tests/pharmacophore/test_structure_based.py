from openpharmacophore import StructureBasedPharmacophore, PharmacophoricPoint
import nglview as nv
import pyunitwizard as puw
import pytest
from copy import deepcopy
import os


def test_init_structured_based_pharmacophore():

    pharmacophore = StructureBasedPharmacophore(None)
    assert pharmacophore.num_frames == 0


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
