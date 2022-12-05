from openpharmacophore import LigandBasedPharmacophore, PharmacophoricPoint
import openpharmacophore.data as data
from matplotlib.colors import to_rgb
import numpy as np
import nglview as nv
import pyunitwizard as puw
import pytest
from rdkit import Chem
from copy import deepcopy
from unittest.mock import Mock, call


def test_init_ligand_based_pharmacophore():
    pharmacophore = LigandBasedPharmacophore()
    assert len(pharmacophore) == 0
    assert len(pharmacophore.ligands) == 0


def test_load_ligands_from_file():
    pharmacophore = LigandBasedPharmacophore()
    pharmacophore.load_ligands(data.ligands["clique_detection.smi"])
    assert len(pharmacophore.ligands) == 5
    assert all([isinstance(lig, Chem.Mol) for lig in pharmacophore.ligands])


def test_is_ligand_file():
    assert LigandBasedPharmacophore._is_ligand_file("mols.smi")
    assert not LigandBasedPharmacophore._is_ligand_file("mols.jpg")


@pytest.fixture()
def pharmacophore_two_ligands():
    pharmacophore = LigandBasedPharmacophore()
    pharmacophore.load_ligands_from_smi([
        "CN[C@H](Cc1ccccc1)C(=O)N2CCC[C@H]2C(=O)NCC3CCC(CC3)N",
        "c1ccc(cc1)S(=O)(=O)CCN2C(=O)N3CC=C[C@H](N3C2=O)C(=O)NC4CCC(CC4)c5cnc([nH]5)N",
    ])
    return pharmacophore


def test_load_ligands_from_smi(pharmacophore_two_ligands):
    assert len(pharmacophore_two_ligands.ligands) == 2
    assert all([isinstance(lig, Chem.Mol) for lig in pharmacophore_two_ligands.ligands])


def test_add_hydrogens_all(pharmacophore_two_ligands):
    pharma = deepcopy(pharmacophore_two_ligands)
    pharma.add_hydrogens("all")
    assert len(
        [a for a in pharma.ligands[0].GetAtoms() if a.GetSymbol() == "H"]) > 0
    assert len(
        [a for a in pharma.ligands[1].GetAtoms() if a.GetSymbol() == "H"]) > 0


def test_add_hydrogens_selection(pharmacophore_two_ligands):
    pharma = deepcopy(pharmacophore_two_ligands)
    pharma.add_hydrogens([1])
    assert len(
        [a for a in pharma.ligands[0].GetAtoms() if a.GetSymbol() == "H"]) == 0
    assert len(
        [a for a in pharma.ligands[1].GetAtoms() if a.GetSymbol() == "H"]) > 0


def test_generate_conformers_all(pharmacophore_two_ligands):
    pharma = deepcopy(pharmacophore_two_ligands)
    pharma.generate_conformers(n_confs=1, ligands="all")
    assert pharma.ligands[0].GetNumConformers() == 1
    assert pharma.ligands[1].GetNumConformers() == 1

    pharma = deepcopy(pharmacophore_two_ligands)
    pharma.generate_conformers(n_confs=[1, 2], ligands="all")
    assert pharma.ligands[0].GetNumConformers() == 1
    assert pharma.ligands[1].GetNumConformers() == 2


def test_generate_conformers_selection(pharmacophore_two_ligands):
    pharma = deepcopy(pharmacophore_two_ligands)
    pharma.generate_conformers(n_confs=1, ligands=[1])
    assert pharma.ligands[0].GetNumConformers() == 0
    assert pharma.ligands[1].GetNumConformers() == 1


def test_generate_conformers_incorrect_selection_len_raises_error(pharmacophore_two_ligands):
    pharma = pharmacophore_two_ligands
    with pytest.raises(ValueError):
        pharma.generate_conformers(n_confs=[1, 2, 3], ligands="all")


def test_find_chem_feats(mocker, pharmacophore_two_ligands):
    mocker.patch(
        "openpharmacophore.pharmacophore.ligand_based.ligand_based.feature_indices",
        side_effect=[
            [(1, 2, 3)],  # ring
            [(4,)],  # hyd
            [(5,)],  # neg charge
            [],  # pos charge
            [(6,), (7,)],  # acceptor
            [],  # donor
            [(1, 2, 3)],  # ring
            [(4,)],  # hyd
            [(5,)],  # neg charge
            [],  # pos charge
            [(6,), (7,)],  # acceptor
            []  # donor
        ]
    )
    pharma = pharmacophore_two_ligands
    pharma.find_chem_feats()
    assert len(pharma.feats) == 2
    expected_feats = {
        "R": [(1, 2, 3)],
        "H": [(4,)],
        "N": [(5,)],
        "P": [],
        "A": [(6,), (7,)],
        "D": [],
    }
    assert pharma.feats[0] == expected_feats
    assert pharma.feats[1] == expected_feats


@pytest.fixture()
def pharmacophore_three_points():
    radius = puw.quantity(1.0, "angstroms")
    center = puw.quantity([1.0, 1.0, 1.0], "angstroms")
    pharma = [
        PharmacophoricPoint("hb donor", center, radius),
        PharmacophoricPoint("aromatic ring", center * 2.0, radius),
        PharmacophoricPoint("hydrophobicity", center * -2.0, radius),
    ]
    pharmacophore = LigandBasedPharmacophore()
    pharmacophore.add_pharmacophore(pharma)
    return pharmacophore


def test_pharmacophore_len(pharmacophore_three_points):
    assert len(pharmacophore_three_points) == 1
    assert len(pharmacophore_three_points[0]) == 3


def test_add_point_to_pharmacophore(pharmacophore_three_points):
    ph = deepcopy(pharmacophore_three_points)
    ph.add_point(PharmacophoricPoint(
        "positive charge",
        puw.quantity([1., 2., 3.], "angstroms"),
        puw.quantity(1., "angstroms")
    ), pharma=0)
    assert len(ph[0]) == 4
    assert ph[0][3].feature_name == "positive charge"


def test_remove_point_from_pharmacophore(pharmacophore_three_points):
    ph = deepcopy(pharmacophore_three_points)
    ph.remove_point(0, 1)
    assert len(ph[0]) == 2
    assert ph[0][0].feature_name == "hb donor"
    assert ph[0][1].feature_name == "hydrophobicity"


def test_pharmacophore_get_item(pharmacophore_three_points):
    assert pharmacophore_three_points[0][0].feature_name == "hb donor"
    assert pharmacophore_three_points[0][1].feature_name == "aromatic ring"


def test_to_rdkit(pharmacophore_three_points):
    rdkit_ph = pharmacophore_three_points.to_rdkit(0)

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
    ), 0)
    ph.add_to_view(view, 0)
    assert len(view._ngl_component_ids) == 5


def setup_pharmacophore_to_file_test(func, file_name, mocker):
    mock_open = mocker.patch(
        "openpharmacophore.pharmacophore.ligand_based.ligand_based.open",
        new=mocker.mock_open())
    func(file_name)
    mock_open.assert_called_once_with(file_name, "w")


def test_to_json(mocker, pharmacophore_three_points):
    file_name = "ph.json"
    setup_pharmacophore_to_file_test(
        pharmacophore_three_points.to_json, file_name, mocker)


def test_to_ligand_scout(mocker, pharmacophore_three_points):
    mock_tree = mocker.patch(
        "openpharmacophore.pharmacophore.ligand_based.ligand_based.io.ligandscout_xml_tree"
    )
    file_name = "ph.pml"
    pharmacophore_three_points.to_ligand_scout(file_name)
    mock_tree.return_value.write.assert_called_once_with(
        file_name, encoding="UTF-8", xml_declaration=True
    )


def test_to_moe(mocker, pharmacophore_three_points):
    file_name = "ph.ph4"
    setup_pharmacophore_to_file_test(
        pharmacophore_three_points.to_moe, file_name, mocker)


def test_to_mol2(mocker, pharmacophore_three_points):
    file_name = "ph.mol2"
    setup_pharmacophore_to_file_test(
        pharmacophore_three_points.to_mol2, file_name, mocker)


def test_pharmacophore_string_representation(pharmacophore_three_points):
    assert pharmacophore_three_points.__repr__() == "LigandBasedPharmacophore(n_pharmacophores: 1)"


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
    pharmacophore.remove_picked_point(mock_view, 0)
    assert len(pharmacophore[0]) == 2
    assert pharmacophore[0][0].short_name == "R"
    assert pharmacophore[0][1].short_name == "H"


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
    pharmacophore.remove_picked_point(mock_view, 0)
    assert len(pharmacophore[0]) == 3


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
        radius=puw.quantity(2.0, "angstroms"),
        pharma=0
    )
    assert len(pharmacophore[0]) == 3
    assert np.all(puw.get_value(pharmacophore[0][0].center) == np.array([1.5, 1.5, 1.5]))
    assert puw.get_value(pharmacophore[0][0].radius) == 2.0


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
        mock_view, "hb acceptor", radius, 0)
    assert len(pharmacophore[0]) == 4
    assert pharmacophore[0][3].short_name == "A"


def test_add_ligands_to_view():

    mock_view = Mock()
    pharmacophore = LigandBasedPharmacophore()
    pharmacophore.load_ligands(data.ligands["clique_detection.smi"])
    pharmacophore.add_ligands_to_view(mock_view)
    ligands_call = [call(lig) for lig in pharmacophore.ligands]

    assert mock_view.add_component.call_count == 5
    assert mock_view.add_component.call_args_list == ligands_call


def ligand_list():
    supplier = Chem.SmilesMolSupplier(data.ligands["clique_detection.smi"],
                                      titleLine=True)
    return [mol for mol in supplier]


def test_show_no_ligands(pharmacophore_three_points):

    view = pharmacophore_three_points.show(ligands=False)
    assert len(view._ngl_component_ids) == 3


def test_show_with_ligands(pharmacophore_three_points):

    pharmacophore = deepcopy(pharmacophore_three_points)
    pharmacophore.ligands = ligand_list()
    view = pharmacophore.show(ligands=True)
    assert len(view._ngl_component_ids) == 8


def test_atom_highlights(pharmacophore_two_ligands):
    pharma = pharmacophore_two_ligands
    pharma._feats = [
        {"R": [(1, 2, 3)]},
        {"D": [(1,)]},
    ]
    atoms, colors, radii = pharma._atom_highlights()

    ring_color = to_rgb(PharmacophoricPoint.palette["aromatic ring"])
    donor_color = to_rgb(PharmacophoricPoint.palette["hb donor"])
    expected_atoms = [(1, 2, 3), (1,)]
    expected_colors = [
        {1: ring_color, 2: ring_color, 3: ring_color},
        {1: donor_color}
    ]
    expected_radii = [
        {1: 0.5, 2: 0.5, 3: 0.5},
        {1: 0.5}
    ]
    assert atoms == expected_atoms
    assert expected_colors == colors
    assert radii == expected_radii


def test_drawing_size():
    assert LigandBasedPharmacophore._drawing_size(300, 280, 2) == (600, 280)
    assert LigandBasedPharmacophore._drawing_size(300, 280, 4) == (1200, 280)
    assert LigandBasedPharmacophore._drawing_size(300, 280, 6) == (1200, 560)
