from openpharmacophore import LigandReceptorPharmacophore, PharmacophoricPoint
import openpharmacophore.data as data
import nglview as nv
import numpy as np
import pyunitwizard as puw
import pytest
from copy import deepcopy
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


# TODO:
#  This mock object is set up so that the distance between the
#  first center of each feature in the ligand is at an optimal distance
#  of the second center of each feature in the protein to create a pharmacophoric
#  point. Right now only distances are being considered. However, angles and other
#  considerations must be made depending on the chemical feature.

ligand_feats_side = [{
        "hb acceptor": [puw.quantity(np.array([0.] * 3), "angstroms"),
                        puw.quantity(np.array([-7.] * 3), "angstroms"),
                        ],
        "aromatic ring": [puw.quantity(np.array([0.] * 3), "angstroms"),
                          puw.quantity(np.array([-7.] * 3), "angstroms"),
                          ],
        "hb donor": [puw.quantity(np.array([0.] * 3), "angstroms"),
                     puw.quantity(np.array([-7.] * 3), "angstroms"),
                     ],
        "hydrophobicity": [puw.quantity(np.array([0.] * 3), "angstroms"),
                           puw.quantity(np.array([-7.] * 3), "angstroms"),
                           ],
        "positive charge": [puw.quantity(np.array([0.] * 3), "angstroms"),
                            puw.quantity(np.array([-7.] * 3), "angstroms"),
                            ],
        "negative charge": [puw.quantity(np.array([0.] * 3), "angstroms"),
                            puw.quantity(np.array([-7.] * 3), "angstroms"),
                            ]
    }]

receptor_feats_side = [{
        "hb acceptor": [puw.quantity(np.array([7.] * 3), "angstroms"),
                        puw.quantity(np.array([2.] * 3), "angstroms"),
                        ],
        "aromatic ring": [puw.quantity(np.array([7.] * 3), "angstroms"),
                          puw.quantity(np.array([2.] * 3), "angstroms"),
                          ],
        "hb donor": [puw.quantity(np.array([7.] * 3), "angstroms"),
                     puw.quantity(np.array([2.] * 3), "angstroms"),
                     ],
        "hydrophobicity": [puw.quantity(np.array([7.] * 3), "angstroms"),
                           puw.quantity(np.array([2.] * 3), "angstroms"),
                           ],
        "positive charge": [puw.quantity(np.array([7.] * 3), "angstroms"),
                            puw.quantity(np.array([2.] * 3), "angstroms"),
                            ],
        "negative charge": [puw.quantity(np.array([7.] * 3), "angstroms"),
                            puw.quantity(np.array([2.] * 3), "angstroms"),
                            ]
    }]


def empty_pharmacophore_one_frame():
    pharmacophore = LigandReceptorPharmacophore()
    pharmacophore.add_frame()
    return pharmacophore


@pytest.mark.skip(reason="Not implemented yet")
def test_hb_acceptor_pharmacophoric_points():
    pharmacophore = empty_pharmacophore_one_frame()
    pharmacophore._hb_acceptor_pharmacophoric_points(
        ligand_feats_side[0]["hb acceptor"],
        receptor_feats_side[0]["hb acceptor"],
        frame=0
    )
    assert len(pharmacophore[0]) == 1
    assert pharmacophore[0][0].feature_name == "hb acceptor"
    assert np.all(puw.get_value(pharmacophore[0][0].center) == np.array([0.]*3))
    assert pharmacophore[0][0].has_direction


@pytest.mark.skip(reason="Not implemented yet")
def test_hb_donor_pharmacophoric_points():
    pharmacophore = empty_pharmacophore_one_frame()
    pharmacophore._hb_donor_pharmacophoric_points(
        ligand_feats_side[0]["hb donor"],
        receptor_feats_side[0]["hb donor"],
        frame=0
    )
    assert len(pharmacophore[0]) == 1
    assert pharmacophore[0].feature_name == "hb donor"
    assert np.all(puw.get_value(pharmacophore[0][0].center) == np.array([0.]*3))
    assert pharmacophore[0][0].has_direction


@pytest.mark.skip(reason="Not implemented yet")
def test_aromatic_pharmacophoric_points():
    pharmacophore = empty_pharmacophore_one_frame()
    pharmacophore._aromatic_pharmacophoric_points(
        ligand_feats_side[0]["aromatic ring"],
        receptor_feats_side[0]["aromatic ring"],
        frame=0
    )
    assert len(pharmacophore[0]) == 1
    assert pharmacophore[0].feature_name == "aromatic ring"
    assert np.all(puw.get_value(pharmacophore[0][0].center) == np.array([0.] * 3))
    assert pharmacophore[0][0].has_direction


def test_hydrophobic_pharmacophoric_points():
    pharmacophore = empty_pharmacophore_one_frame()
    pharmacophore._hydrophobic_pharmacophoric_points(
        ligand_feats_side[0]["hydrophobicity"],
        receptor_feats_side[0]["hydrophobicity"],
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
    assert np.allclose(puw.get_value(points[0].center), np.array([1/3, 2/3, 0]))
    assert np.allclose(puw.get_value(points[1].center), np.array([1/3, 3/2, 0]))


def test_charge_pharmacophoric_points():
    pharmacophore = empty_pharmacophore_one_frame()
    pharmacophore._charge_pharmacophoric_points(
        ligand_feats_side[0]["positive charge"],
        receptor_feats_side[0]["negative charge"],
        "positive charge",
        frame=0
    )
    assert len(pharmacophore[0]) == 1
    assert pharmacophore[0][0].feature_name == "positive charge"
    assert np.all(puw.get_value(pharmacophore[0][0].center) == np.array([0.] * 3))
    assert not pharmacophore[0][0].has_direction


@pytest.mark.skip(reason="Not implemented yet")
def test_extract_all_features(mocker):
    pharmacophore = LigandReceptorPharmacophore()
    pharmacophore._pl_complex = mocker.Mock()
    pharmacophore._pl_complex.ligand_feats_center.side_effect = ligand_feats_side
    pharmacophore._pl_complex.receptor_feats_center.side_effect = receptor_feats_side

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
