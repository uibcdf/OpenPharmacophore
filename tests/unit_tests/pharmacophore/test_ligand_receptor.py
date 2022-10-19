from openpharmacophore import LigandReceptorPharmacophore, PharmacophoricPoint
import openpharmacophore.data as data
import nglview as nv
import numpy as np
import pyunitwizard as puw
import pytest
from copy import deepcopy
from collections import namedtuple
import os


def test_init_ligand_receptor_pharmacophore():
    pharmacophore = LigandReceptorPharmacophore()
    assert len(pharmacophore) == 0


def test_load_pdb_file(mocker):
    mock_pl_complex = mocker.patch(
        "openpharmacophore.pharmacophore.ligand_receptor.PLComplex")
    pharmacophore = LigandReceptorPharmacophore()
    pharmacophore.load_pdb(data.pdb["1ncr.pdb"])
    assert pharmacophore.num_frames == 1
    assert len(pharmacophore) == 0

    mock_pl_complex.assert_called_once_with(data.pdb["1ncr.pdb"])


def test_load_pdb_id(mocker):
    mocker.patch("openpharmacophore.LigandReceptorPharmacophore._fetch_pdb",
                 return_value=b"pdb")
    mock_pl_complex = mocker.patch(
        "openpharmacophore.pharmacophore.ligand_receptor.PLComplex")
    pharmacophore = LigandReceptorPharmacophore()
    pharmacophore.load_pdb_id("1NCR")
    assert len(pharmacophore) == 0
    assert pharmacophore.num_frames == 1
    # TODO: check that PLComplex is instantiated with the name of the temporary file
    mock_pl_complex.assert_called_once()


def test_is_pdb_id():
    assert LigandReceptorPharmacophore._is_pdb_id("1A52")
    assert LigandReceptorPharmacophore._is_pdb_id("5UFW")
    assert not LigandReceptorPharmacophore._is_pdb_id("1A52D")
    assert not LigandReceptorPharmacophore._is_pdb_id("1A3")


def test_fetch_pdb(mocker):
    mock_get = mocker.patch("requests.get")
    mock_get.return_value.status_code = 200
    mock_get.return_value.content = b"pdb"

    assert LigandReceptorPharmacophore._fetch_pdb("1A52") == b"pdb"
    mock_get.assert_called_once_with('http://files.rcsb.org/download/1A52.pdb', allow_redirects=True)


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

    ph = LigandReceptorPharmacophore()
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
    rdkit_ph = pharmacophore_one_frame.to_rdkit(0)
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


ChemFeats = namedtuple("ChemFeats", [
    "acc_cent", "acc_ind", "don_cent", "don_ind",
    "aro_cent", "aro_ind", "hyd_cent", "charge_cent",
])


@pytest.fixture()
def ligand_chem_feats():
    acceptors_centers = [puw.quantity(np.array([0.] * 3), "angstroms"),
                         puw.quantity(np.array([-7.] * 3), "angstroms")]
    acceptor_indices = [[1], [1]]

    donors_centers = [puw.quantity(np.array([0.] * 3), "angstroms"),
                      puw.quantity(np.array([-7.] * 3), "angstroms")]
    donors_indices = [[1], [1]]

    aromatic_centers = [puw.quantity(np.array([0., 0., 0.]), "angstroms")]
    aromatic_indices = [[0, 1, 2, 3, 4, 5]]

    hyd_centers = [puw.quantity(np.array([0.] * 3), "angstroms"),
                   puw.quantity(np.array([-7.] * 3), "angstroms")]

    charge_centers = [puw.quantity(np.array([0.] * 3), "angstroms"),
                      puw.quantity(np.array([-7.] * 3), "angstroms")]

    return ChemFeats(acceptors_centers, acceptor_indices,
                     donors_centers, donors_indices,
                     aromatic_centers, aromatic_indices,
                     hyd_centers, charge_centers)


@pytest.fixture()
def receptor_chem_feats():
    acceptors_centers = [puw.quantity(np.array([7.] * 3), "angstroms"),
                         puw.quantity(np.array([2.] * 3), "angstroms")]
    acceptor_indices = [[1], [1]]

    aromatic_centers = [puw.quantity(np.array([1., 0., 3.]), "angstroms")]
    aromatic_indices = [[6, 7, 8, 9, 10, 11]]

    donors_centers = [puw.quantity(np.array([7.] * 3), "angstroms"),
                      puw.quantity(np.array([2.] * 3), "angstroms")]
    donors_indices = [[1], [1]]

    hyd_centers = [puw.quantity(np.array([7.] * 3), "angstroms"),
                   puw.quantity(np.array([2.] * 3), "angstroms")]

    charge_centers = [puw.quantity(np.array([7.] * 3), "angstroms"),
                      puw.quantity(np.array([2.] * 3), "angstroms")]

    return ChemFeats(acceptors_centers, acceptor_indices,
                     donors_centers, donors_indices,
                     aromatic_centers, aromatic_indices,
                     hyd_centers, charge_centers)


def empty_pharmacophore_one_frame():
    pharmacophore = LigandReceptorPharmacophore()
    pharmacophore.add_frame()
    return pharmacophore


@pytest.mark.skip(reason="Not implemented yet")
def test_hb_acceptor_pharmacophoric_points(ligand_chem_feats, receptor_chem_feats):
    pharmacophore = empty_pharmacophore_one_frame()
    pharmacophore._hb_acceptor_pharmacophoric_points(
        ligand_chem_feats.acc_cent,
        ligand_chem_feats.acc_ind,
        receptor_chem_feats.don_cent,
        receptor_chem_feats.don_ind,
        frame=0
    )
    assert len(pharmacophore[0]) == 1
    assert pharmacophore[0][0].feature_name == "hb acceptor"
    assert np.all(puw.get_value(pharmacophore[0][0].center) == np.array([0.] * 3))
    assert pharmacophore[0][0].has_direction


@pytest.mark.skip(reason="Not implemented yet")
def test_hb_donor_pharmacophoric_points(ligand_chem_feats, receptor_chem_feats):
    pharmacophore = empty_pharmacophore_one_frame()
    pharmacophore._hb_donor_pharmacophoric_points(
        ligand_chem_feats.don_center,
        ligand_chem_feats.don_ind,
        receptor_chem_feats.acc_cent,
        receptor_chem_feats.acc_ind,
        frame=0
    )
    assert len(pharmacophore[0]) == 1
    assert pharmacophore[0].feature_name == "hb donor"
    assert np.all(puw.get_value(pharmacophore[0][0].center) == np.array([0.] * 3))
    assert pharmacophore[0][0].has_direction


def test_aromatic_pharmacophoric_points_exceeds_max_distance():
    pharmacophore = empty_pharmacophore_one_frame()
    lig_centers = [puw.quantity(np.array([0., 0., 0.]), "angstroms")]
    lig_indices = [[6, 7, 8, 9, 10, 11]]
    rec_centers = [puw.quantity(np.array([10., 10., 10.]), "angstroms")]
    rec_indices = [0, 1, 2, 3, 4, 5]
    pharmacophore._aromatic_pharmacophoric_points(
        lig_centers, lig_indices,
        rec_centers, rec_indices,
        frame=0
    )
    assert len(pharmacophore[0]) == 0


def test_aromatic_pharmacophoric_points_pstack(mocker, ligand_chem_feats, receptor_chem_feats):
    pharmacophore = empty_pharmacophore_one_frame()
    pharmacophore._pl_complex = mocker.Mock()
    pharmacophore._pl_complex.coords = np.array([
        # First ring is a hexagon on the xy plane
        [1, 0, 0],
        [1/2, np.sqrt(3)/2, 0],
        [-1/2, np.sqrt(3)/2, 0],
        [-1, 0, 0],
        [-1/2, -np.sqrt(3)/2, 0],
        [1/2, -np.sqrt(3)/2, 0],
        # Second ring is a hexagon on the plane z = 3
        [2, 0, 3],
        [3/2, np.sqrt(3)/2, 3],
        [1/2, np.sqrt(3)/2, 3],
        [0, 0, 3],
        [1/2, -np.sqrt(3)/2, 3],
        [3/2, -np.sqrt(3)/2, 3],
    ])
    pharmacophore._aromatic_pharmacophoric_points(
        ligand_chem_feats.aro_cent,
        ligand_chem_feats.aro_ind,
        receptor_chem_feats.aro_cent,
        receptor_chem_feats.aro_ind,
        frame=0
    )
    assert len(pharmacophore[0]) == 1
    assert pharmacophore[0][0].feature_name == "aromatic ring"
    assert np.all(puw.get_value(pharmacophore[0][0].center) == np.array([0.] * 3))
    assert pharmacophore[0][0].has_direction
    assert np.allclose(pharmacophore[0].direction,
                       np.array([1/np.sqrt(10), 0, 3/np.sqrt(10)]))


def test_aromatic_pharmacophoric_points_exceeds_max_offset(mocker, ligand_chem_feats, receptor_chem_feats):
    pharmacophore = empty_pharmacophore_one_frame()
    pharmacophore._pl_complex = mocker.Mock()
    pharmacophore._pl_complex.coords = np.array([
        # First ring is a hexagon on the xy plane
        [1, 0, 0],
        [1 / 2, np.sqrt(3) / 2, 0],
        [-1 / 2, np.sqrt(3) / 2, 0],
        [-1, 0, 0],
        [-1 / 2, -np.sqrt(3) / 2, 0],
        [1 / 2, -np.sqrt(3) / 2, 0],
        # Second ring is a hexagon on the plane z = 3
        [6, 0, 3],
        [11 / 2, np.sqrt(3) / 2, 3],
        [9 / 2, np.sqrt(3) / 2, 3],
        [4, 0, 3],
        [11 / 2, -np.sqrt(3) / 2, 3],
        [9 / 2, -np.sqrt(3) / 2, 3],
    ])
    pharmacophore._aromatic_pharmacophoric_points(
        ligand_chem_feats.aro_cent,
        ligand_chem_feats.aro_ind,
        [puw.quantity(np.array([5, 0, 3]), "angstroms")],
        receptor_chem_feats.aro_ind,
        frame=0
    )
    assert len(pharmacophore[0]) == 0


def test_aromatic_pharmacophoric_points_t_stack(mocker, ligand_chem_feats, receptor_chem_feats):
    pharmacophore = empty_pharmacophore_one_frame()
    pharmacophore._pl_complex = mocker.Mock()
    pharmacophore._pl_complex.coords = np.array([
        # First ring is a hexagon on the xy plane
        [1, 0, 0],
        [1 / 2, np.sqrt(3) / 2, 0],
        [-1 / 2, np.sqrt(3) / 2, 0],
        [-1, 0, 0],
        [-1 / 2, -np.sqrt(3) / 2, 0],
        [1 / 2, -np.sqrt(3) / 2, 0],
        # Second ring is a hexagon on the plane x = 2
        [2, 1, 0],
        [2, 1/2, np.sqrt(3) / 2],
        [2, -1/2, np.sqrt(3) / 2],
        [2, -1, 0],
        [2, -1/2, -np.sqrt(3) / 2],
        [2, 1/2, -np.sqrt(3) / 2],
    ])
    pharmacophore._aromatic_pharmacophoric_points(
        ligand_chem_feats.aro_cent,
        ligand_chem_feats.aro_ind,
        [puw.quantity(np.array([2, 0, 0]), "angstroms")],
        receptor_chem_feats.aro_ind,
        frame=0
    )
    assert len(pharmacophore[0]) == 1
    assert pharmacophore[0][0].feature_name == "aromatic ring"
    assert np.all(puw.get_value(pharmacophore[0][0].center) == np.array([0.] * 3))
    assert pharmacophore[0][0].has_direction
    assert np.allclose(pharmacophore[0].direction,
                       np.array([0., 0., 1.]))


def test_hydrophobic_pharmacophoric_points(ligand_chem_feats, receptor_chem_feats):
    pharmacophore = empty_pharmacophore_one_frame()
    pharmacophore._hydrophobic_pharmacophoric_points(
        ligand_chem_feats.hyd_cent,
        receptor_chem_feats.hyd_cent,
        frame=0
    )
    assert len(pharmacophore[0]) == 1
    assert pharmacophore[0][0].feature_name == "hydrophobicity"
    assert np.all(puw.get_value(pharmacophore[0][0].center) == np.array([0.] * 3))
    assert not pharmacophore[0][0].has_direction


def test_merge_hydrophobic_points():
    radius = puw.quantity(1.0, "angstroms")
    centers = [
        puw.quantity(np.array([0., 0., 0.]), "angstroms"),
        puw.quantity(np.array([1., 1., 0.]), "angstroms"),
        puw.quantity(np.array([0., 1., 0.]), "angstroms"),
        puw.quantity(np.array([0., 2.5, 0.]), "angstroms"),
    ]
    points = LigandReceptorPharmacophore._merge_hydrophobic_points(
        centers, radius
    )
    assert len(points) == 2
    assert all([p.feature_name == "hydrophobicity" for p in points])
    assert np.allclose(puw.get_value(points[0].center), np.array([1 / 3, 2 / 3, 0]))
    assert np.allclose(puw.get_value(points[1].center), np.array([1 / 3, 3 / 2, 0]))


def test_charge_pharmacophoric_points(ligand_chem_feats, receptor_chem_feats):
    pharmacophore = empty_pharmacophore_one_frame()
    pharmacophore._charge_pharmacophoric_points(
        ligand_chem_feats.charge_cent,
        receptor_chem_feats.charge_cent,
        "positive charge",
        frame=0
    )
    assert len(pharmacophore[0]) == 1
    assert pharmacophore[0][0].feature_name == "positive charge"
    assert np.all(puw.get_value(pharmacophore[0][0].center) == np.array([0.] * 3))
    assert not pharmacophore[0][0].has_direction


@pytest.mark.skip(reason="Not implemented yet")
def test_extract_all_features(mocker, ligand_chem_feats, receptor_chem_feats):
    pharmacophore = LigandReceptorPharmacophore()
    pharmacophore._pl_complex = mocker.Mock()

    pharmacophore.extract("EST:B")
    assert len(pharmacophore[0]) == 6
    expected_features = {
        "hb acceptor",
        "aromatic ring",
        "hb donor",
        "hydrophobicity",
        "positive charge",
        "negative charge",
    }
    actual_features = set()
    for pnt in pharmacophore[0]:
        actual_features.add(pnt.feature_name)
    assert actual_features == expected_features

    center_expected = np.array([0.] * 3)
    for point in pharmacophore[0]:
        assert np.all(point.center == center_expected)
