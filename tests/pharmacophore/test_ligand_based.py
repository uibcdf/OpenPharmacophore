from openpharmacophore import LigandBasedPharmacophore, PharmacophoricPoint
import openpharmacophore.data as data
import numpy as np
import nglview as nv
import pyunitwizard as puw
import pytest
from rdkit import Chem
from copy import deepcopy
import os
from unittest.mock import Mock, call


def test_init_with_pharmacophore_file():
    pharmacophore = LigandBasedPharmacophore(data.pharmacophores["ligscout"])
    assert len(pharmacophore) == 4

    pharmacophore = LigandBasedPharmacophore(data.pharmacophores["gmp"])
    assert len(pharmacophore) == 10


def test_init_with_smi_file():
    ligands = data.ligands["clique_detection"]
    pharmacophore = LigandBasedPharmacophore(ligands)
    assert len(pharmacophore) == 0
    assert len(pharmacophore.ligands) == 5


def test_init_with_mol2_file():
    ligands = data.ligands["ace"]
    pharmacophore = LigandBasedPharmacophore(ligands)
    assert len(pharmacophore) == 0
    assert len(pharmacophore.ligands) == 3


def test_init_with_sdf_file():
    ligands = data.ligands["sdf_example"]
    pharmacophore = LigandBasedPharmacophore(ligands)
    assert len(pharmacophore) == 0
    assert len(pharmacophore.ligands) == 3


def test_from_file_moe():
    pharmacophore = LigandBasedPharmacophore([])
    pharmacophore.from_file(data.pharmacophores["gmp"])
    assert len(pharmacophore) == 10


def test_from_file_mol2():
    pharmacophore = LigandBasedPharmacophore([])
    pharmacophore.from_file(data.pharmacophores["elastase"])
    assert len(pharmacophore) == 4


def test_from_file_json():
    pharmacophore = LigandBasedPharmacophore([])
    pharmacophore.from_file(data.pharmacophores["1M70"])
    assert len(pharmacophore) == 5


def test_is_ligand_file():
    assert LigandBasedPharmacophore._is_ligand_file("mols.smi")
    assert not LigandBasedPharmacophore._is_ligand_file("mols.jpg")


@pytest.fixture()
def pharmacophore_three_points():
    radius = puw.quantity(1.0, "angstroms")
    center = puw.quantity([1.0, 1.0, 1.0], "angstroms")
    points = [
        PharmacophoricPoint("hb donor", center, radius),
        PharmacophoricPoint("aromatic ring", center * 2.0, radius),
        PharmacophoricPoint("hydrophobicity", center * -2.0, radius),
    ]
    pharmacophore = LigandBasedPharmacophore([])
    pharmacophore.pharmacophoric_points = points
    return pharmacophore


def test_pharmacophore_len(pharmacophore_three_points):
    assert len(pharmacophore_three_points) == 3


def test_iterate_pharmacophore(pharmacophore_three_points):
    short_names = ["D", "R", "H"]
    ii = 0
    iterated = False
    for point in pharmacophore_three_points:
        assert point.short_name == short_names[ii]
        ii += 1
        iterated = True

    assert iterated


def test_add_point_to_pharmacophore(pharmacophore_three_points):
    ph = deepcopy(pharmacophore_three_points)
    ph.add_point(PharmacophoricPoint(
        "positive charge",
        puw.quantity([1., 2., 3.], "angstroms"),
        puw.quantity(1., "angstroms")
    ))
    assert len(ph) == 4
    assert ph[3].feature_name == "positive charge"


def test_remove_point_from_pharmacophore(pharmacophore_three_points):
    ph = deepcopy(pharmacophore_three_points)
    ph.remove_point(1)
    assert len(ph) == 2
    assert ph[0].feature_name == "hb donor"
    assert ph[1].feature_name == "hydrophobicity"


def test_pharmacophore_get_item(pharmacophore_three_points):
    assert pharmacophore_three_points[0].feature_name == "hb donor"
    assert pharmacophore_three_points[1].feature_name == "aromatic ring"
    assert pharmacophore_three_points[2].feature_name == "hydrophobicity"


def test_pharmacophore_equality(pharmacophore_three_points):
    radius = puw.quantity(1.0, "angstroms")
    donor = PharmacophoricPoint("hb donor",
                                puw.quantity([1.0, 1.0, 1.0], "angstroms"),
                                radius)
    pharmacophore_1 = LigandBasedPharmacophore([])
    pharmacophore_1.pharmacophoric_points = [donor]
    assert not pharmacophore_1 == pharmacophore_three_points

    ring = PharmacophoricPoint("aromatic ring",
                               puw.quantity([2.0, 2.0, 2.0], "angstroms"),
                               radius)
    hydrophobic = PharmacophoricPoint("hydrophobicity",
                                      puw.quantity([-1.0, 2.0, 2.0], "angstroms"),
                                      radius)
    pharmacophore_2 = LigandBasedPharmacophore([])
    pharmacophore_2.pharmacophoric_points = [donor, ring, hydrophobic]
    assert not pharmacophore_2 == pharmacophore_three_points

    hydrophobic.center = puw.quantity([-2.0, -2.0, -2.0], "angstroms")
    pharmacophore_3 = LigandBasedPharmacophore([])
    pharmacophore_3.pharmacophoric_points = [donor, ring, hydrophobic]
    assert pharmacophore_3 == pharmacophore_three_points


def test_to_rdkit(pharmacophore_three_points):
    rdkit_ph, radii = pharmacophore_three_points.to_rdkit()

    assert len(radii) == 3
    assert np.allclose(np.array(radii), np.array([1.0, 1.0, 1.0]))

    feats = rdkit_ph.getFeatures()
    assert len(feats) == 3

    acceptor = feats[0]
    assert acceptor.GetFamily() == "Donor"
    assert np.allclose(acceptor.GetPos().x, 1.0)
    assert np.allclose(acceptor.GetPos().y, 1.0)
    assert np.allclose(acceptor.GetPos().z, 1.0)

    ring_1 = feats[1]
    assert ring_1.GetFamily() == "Aromatic"
    assert np.allclose(ring_1.GetPos().x, 2.0)
    assert np.allclose(ring_1.GetPos().y, 2.0)
    assert np.allclose(ring_1.GetPos().z, 2.0)

    ring_2 = feats[2]
    assert ring_2.GetFamily() == "Hydrophobe"
    assert np.allclose(ring_2.GetPos().x, -2.0)
    assert np.allclose(ring_2.GetPos().y, -2.0)
    assert np.allclose(ring_2.GetPos().z, -2.0)


def test_distance_matrix():
    # Test pharmacophore with zero elements
    pharmacophore = LigandBasedPharmacophore([])
    matrix = pharmacophore.distance_matrix()
    assert matrix.shape == (0, 0)

    radius = puw.quantity(1.0, "angstroms")
    points = [
        PharmacophoricPoint(feat_type="hb acceptor",
                            center=puw.quantity([1, 0, -5], "angstroms"),
                            radius=radius),
        PharmacophoricPoint(feat_type="hb donor",
                            center=puw.quantity([2, 1, 0], "angstroms"),
                            radius=radius),
        PharmacophoricPoint(feat_type="aromatic ring",
                            center=puw.quantity([-3, 2, -1], "angstroms"),
                            radius=radius),
    ]

    pharmacophore = LigandBasedPharmacophore([])
    pharmacophore.pharmacophoric_points = points
    distance_matrix = pharmacophore.distance_matrix()

    assert distance_matrix.shape == (3, 3)
    assert np.allclose(distance_matrix,
                       np.array([[0, np.sqrt(27), 6],
                                 [np.sqrt(27), 0, np.sqrt(27)],
                                 [6, np.sqrt(27), 0]])
                       )

    points = [
        PharmacophoricPoint(feat_type="hb acceptor",
                            center=puw.quantity([1, 2, 3], "angstroms"),
                            radius=radius),
        PharmacophoricPoint(feat_type="negative charge",
                            center=puw.quantity([0, 0, 0], "angstroms"),
                            radius=radius),
        PharmacophoricPoint(feat_type="positive charge",
                            center=puw.quantity([-1, 0, -1], "angstroms"),
                            radius=radius),
        PharmacophoricPoint(feat_type="aromatic ring",
                            center=puw.quantity([2, -1, 1], "angstroms"),
                            radius=radius),
    ]

    sq = np.sqrt
    pharmacophore = LigandBasedPharmacophore([])
    pharmacophore.pharmacophoric_points = points
    distance_matrix = pharmacophore.distance_matrix()
    assert distance_matrix.shape == (4, 4)
    assert np.allclose(distance_matrix,
                       np.array([[0, sq(14), sq(24), sq(14)],
                                 [sq(14), 0, sq(2), sq(6)],
                                 [sq(24), sq(2), 0, sq(14)],
                                 [sq(14), sq(6), sq(14), 0]]))


def test_adding_pharmacophore_to_view_updates_components(pharmacophore_three_points):
    view = nv.NGLWidget()
    assert len(view._ngl_component_ids) == 0

    # The view should have a component for each sphere and one for each vector
    ph = deepcopy(pharmacophore_three_points)
    ph.add_point(PharmacophoricPoint(
        "negative charge",
        puw.quantity([2., 3., 1.], "angstroms"),
        puw.quantity(2., "angstroms"),
        direction=[0.3, 0.3, 0.3]
    ))
    ph.add_to_view(view)
    assert len(view._ngl_component_ids) == 5


def assert_file_is_created(file_name):
    assert os.path.isfile(file_name)
    os.remove(file_name)


def test_to_json(pharmacophore_three_points):
    file_name = "ph.json"
    pharmacophore_three_points.to_json(file_name)
    assert_file_is_created(file_name)


def test_to_ligand_scout(pharmacophore_three_points):
    file_name = "ph.pml"
    pharmacophore_three_points.to_ligand_scout(file_name)
    assert_file_is_created(file_name)


def test_to_moe(pharmacophore_three_points):
    file_name = "ph.ph4"
    pharmacophore_three_points.to_moe(file_name)
    assert_file_is_created(file_name)


def test_to_mol2(pharmacophore_three_points):
    file_name = "ph.mol2"
    pharmacophore_three_points.to_mol2(file_name)
    assert_file_is_created(file_name)


def test_pharmacophore_string_representation(pharmacophore_three_points):
    assert pharmacophore_three_points.__repr__() == "LigandBasedPharmacophore(n_pharmacophoric_points: 3)"


def test_remove_picked_point(pharmacophore_three_points):
    mock_view = Mock()
    mock_view._ngl_component_names = ["nglview.adaptor.RdkitStructure",
                                      'nglview.shape.Shape',
                                      'nglview.shape.Shape',
                                      'nglview.shape.Shape']
    # We want to remove the first element from the pharmacophore
    mock_view.picked = {
        "component": 1
    }
    pharmacophore = deepcopy(pharmacophore_three_points)
    pharmacophore.remove_picked_point(mock_view)
    assert len(pharmacophore) == 2
    assert pharmacophore[0].short_name == "R"
    assert pharmacophore[1].short_name == "H"


def test_remove_picked_point_with_no_selected_point(pharmacophore_three_points):
    # Suppose we select an atom from the molecule
    # The pharmacophore should be intact
    mock_view = Mock()
    mock_view._ngl_component_names = ["nglview.adaptor.RdkitStructure",
                                      'nglview.shape.Shape']
    mock_view.picked = {
        "atom1": {
            "x": 1,
            "y": 1,
            "z": 1,
        },
        "component": 0
    }
    pharmacophore = deepcopy(pharmacophore_three_points)
    pharmacophore.remove_picked_point(mock_view)
    assert len(pharmacophore) == 3


def test_edit_picked_point(pharmacophore_three_points):
    mock_view = Mock()
    mock_view._ngl_component_names = ["nglview.adaptor.RdkitStructure",
                                      'nglview.shape.Shape',
                                      'nglview.shape.Shape',
                                      'nglview.shape.Shape']
    # We want to remove the first element from the pharmacophore
    mock_view.picked = {
        "component": 1
    }
    pharmacophore = deepcopy(pharmacophore_three_points)
    pharmacophore.edit_picked_point(
        mock_view,
        center=puw.quantity([1.5, 1.5, 1.5], "angstroms"),
        radius=puw.quantity(2.0, "angstroms"))
    assert len(pharmacophore) == 3
    assert np.all(puw.get_value(pharmacophore[0].center) == np.array([1.5, 1.5, 1.5]))
    assert puw.get_value(pharmacophore[0].radius) == 2.0


def test_add_point_in_picked_location(pharmacophore_three_points):
    mock_view = Mock()
    mock_view._ngl_component_names = ["nglview.adaptor.RdkitStructure",
                                      'nglview.shape.Shape',
                                      'nglview.shape.Shape',
                                      'nglview.shape.Shape']
    mock_view.picked = {
        "atom1": {
            "x": 3.5,
            "y": 3.5,
            "z": 3.5,
        },
        "component": 0
    }

    pharmacophore = deepcopy(pharmacophore_three_points)
    radius = puw.quantity(2.0, "angstroms")
    pharmacophore.add_point_in_picked_location(
        mock_view, "hb acceptor", radius)
    assert len(pharmacophore) == 4
    assert pharmacophore[3].short_name == "A"


def test_load_ligands_smi_file():
    ligands = LigandBasedPharmacophore._load_ligands(
        data.ligands["clique_detection"]
    )
    assert len(ligands) == 5


def test_load_ligands_mol2_file():
    ligands = LigandBasedPharmacophore._load_ligands(
        data.ligands["ace"]
    )
    assert len(ligands) == 3


def test_load_ligands_sfd_file():
    ligands = LigandBasedPharmacophore._load_ligands(
        data.ligands["sdf_example"]
    )
    assert len(ligands) == 3


def test_add_ligands_to_view():

    mock_view = Mock()
    pharmacophore = LigandBasedPharmacophore(data.ligands["clique_detection"])
    pharmacophore.add_ligands_to_view(mock_view)
    ligands_call = [call(lig) for lig in pharmacophore.ligands]

    assert mock_view.add_component.call_count == 5
    assert mock_view.add_component.call_args_list == ligands_call


def ligand_list():
    supplier = Chem.SmilesMolSupplier(data.ligands["clique_detection"],
                                      titleLine=False)
    return [mol for mol in supplier]


def test_show_no_ligands(pharmacophore_three_points):

    view = pharmacophore_three_points.show(ligands=False)
    assert len(view._ngl_component_ids) == 3


def test_show_with_ligands(pharmacophore_three_points):

    pharmacophore = deepcopy(pharmacophore_three_points)
    pharmacophore.ligands = ligand_list()
    view = pharmacophore.show(ligands=True)
    assert len(view._ngl_component_ids) == 8
