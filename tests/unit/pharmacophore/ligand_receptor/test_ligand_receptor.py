from openpharmacophore import LigandReceptorPharmacophore, PharmacophoricPoint, Pharmacophore
import nglview as nv
import numpy as np
import pyunitwizard as puw
import pytest
from copy import deepcopy
from collections import namedtuple

# Imports for mocking
lr_module = "openpharmacophore.pharmacophore.ligand_receptor.ligand_receptor"
lr_class = lr_module + ".LigandReceptorPharmacophore"


def test_init_ligand_receptor_pharmacophore():
    pharmacophore = LigandReceptorPharmacophore()
    assert len(pharmacophore) == 0


def test_load_pdb_file(mocker, pdb_1ncr_path):
    mock_pl_complex = mocker.patch(
        lr_module + ".PLComplex")
    pharmacophore = LigandReceptorPharmacophore()
    pharmacophore.load_receptor(pdb_1ncr_path)
    assert len(pharmacophore) == 0

    mock_pl_complex.assert_called_once_with(pdb_1ncr_path)


def test_load_pdb_id(mocker):
    mocker.patch(lr_class + "._fetch_pdb",
                 return_value=b"pdb")
    mock_pl_complex = mocker.patch(
        lr_module + ".PLComplex")
    pharmacophore = LigandReceptorPharmacophore()
    pharmacophore.load_pdb_id("1NCR")
    assert len(pharmacophore) == 0
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
    ph.add_pharmacophore(
        Pharmacophore(points=[acceptor, donor, aromatic])
    )
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


def setup_pharmacophore_to_file_test(func, file_name, mocker):
    mock_open = mocker.patch(
        "openpharmacophore.pharmacophore.ligand_receptor.ligand_receptor.open",
        new=mocker.mock_open())
    func(file_name, frame=0)
    mock_open.assert_called_once_with(file_name, "w")


def test_to_json(mocker, pharmacophore_one_frame):
    file_name = "ph.json"
    setup_pharmacophore_to_file_test(
        pharmacophore_one_frame.to_json, file_name, mocker)


def test_to_ligand_scout(mocker, pharmacophore_one_frame):
    mock_tree = mocker.patch(
        "openpharmacophore.pharmacophore.ligand_based.ligand_based.io.ligandscout_xml_tree"
    )
    file_name = "ph.pml"
    pharmacophore_one_frame.to_ligand_scout(file_name, frame=0)
    mock_tree.return_value.write.assert_called_once_with(
        file_name, encoding="UTF-8", xml_declaration=True
    )


def test_to_moe(mocker, pharmacophore_one_frame):
    file_name = "ph.ph4"
    setup_pharmacophore_to_file_test(
        pharmacophore_one_frame.to_moe, file_name, mocker)


def test_to_mol2(mocker, pharmacophore_one_frame):
    file_name = "ph.mol2"
    setup_pharmacophore_to_file_test(
        pharmacophore_one_frame.to_mol2, file_name, mocker)


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
    "aro_cent", "aro_ind", "hyd_cent", "charge_cent",
])


@pytest.fixture()
def ligand_chem_feats():
    aromatic_centers = [puw.quantity(np.array([0., 0., 0.]), "angstroms")]
    aromatic_indices = [[0, 1, 2, 3, 4, 5]]

    hyd_centers = [puw.quantity(np.array([0.] * 3), "angstroms"),
                   puw.quantity(np.array([-7.] * 3), "angstroms")]

    charge_centers = [puw.quantity(np.array([0.] * 3), "angstroms"),
                      puw.quantity(np.array([-7.] * 3), "angstroms")]

    return ChemFeats(aromatic_centers, aromatic_indices,
                     hyd_centers, charge_centers)


@pytest.fixture()
def receptor_chem_feats():
    aromatic_centers = [puw.quantity(np.array([1., 0., 3.]), "angstroms")]
    aromatic_indices = [[6, 7, 8, 9, 10, 11]]

    hyd_centers = [puw.quantity(np.array([7.] * 3), "angstroms"),
                   puw.quantity(np.array([2.] * 3), "angstroms")]

    charge_centers = [puw.quantity(np.array([7.] * 3), "angstroms"),
                      puw.quantity(np.array([2.] * 3), "angstroms")]

    return ChemFeats(aromatic_centers, aromatic_indices,
                     hyd_centers, charge_centers)


def empty_pharmacophore_one_frame():
    pharmacophore = LigandReceptorPharmacophore()
    pharmacophore.add_frame()
    return pharmacophore


def set_up_hbond_pharmacophoric_points(mocker):
    pharma = empty_pharmacophore_one_frame()
    pharma._pl_complex = mocker.Mock()
    pharma.receptor.coords = puw.quantity(
        np.array([[float(ii)] * 3 for ii in range(0, 9)]), "angstroms")
    pharma.receptor.coords = pharma.receptor.coords.reshape((1, 9, 3))
    assert pharma.receptor.coords.shape == (1, 9, 3)
    pharma.receptor.lig_indices = [5, 6, 7, 8]

    h_bonds = np.array([
        [0, 2, 5],
        [1, 3, 6],
        [7, 8, 4],
    ])

    return pharma, h_bonds


def test_hbond_donor_pharmacophoric_points(mocker):
    pharma, h_bonds = set_up_hbond_pharmacophoric_points(mocker)
    pharma._hbond_donor_pharmacophoric_points(h_bonds, 0)
    assert len(pharma[0]) == 1
    assert pharma[0][0].feature_name == "hb donor"
    assert np.all(puw.get_value(pharma[0][0].center) ==
                  np.array([7., 7., 7.]))

    # All points lie on the same line, they all have the same direction
    expected_direction = np.array([1 / np.sqrt(3)] * 3) * -1
    assert pharma[0][0].has_direction
    assert np.allclose(pharma[0][0].direction, expected_direction)


def test_hb_acceptor_pharmacophoric_points(mocker):
    pharma, h_bonds = set_up_hbond_pharmacophoric_points(mocker)
    pharma._hbond_acceptor_pharmacophoric_points(h_bonds, 0)
    assert len(pharma[0]) == 2
    assert pharma[0][0].feature_name == "hb acceptor"
    assert np.all(puw.get_value(pharma[0][0].center) ==
                  np.array([5., 5., 5.]))

    assert pharma[0][1].feature_name == "hb acceptor"
    assert np.all(puw.get_value(pharma[0][1].center) ==
                  np.array([6., 6., 6.]))

    expected_direction = np.array([1 / np.sqrt(3)] * 3) * -1
    assert all(p.has_direction for p in pharma[0])
    assert all(np.allclose(pharma[0][ii].direction, expected_direction)
               for ii in range(len(pharma[0])))


def test_aromatic_pharmacophoric_points_exceeds_max_distance():
    pharmacophore = empty_pharmacophore_one_frame()
    lig_centers = [puw.quantity(np.array([0., 0., 0.]), "angstroms")]
    lig_indices = [[6, 7, 8, 9, 10, 11]]
    rec_centers = [puw.quantity(np.array([10., 10., 10.]), "angstroms")]
    rec_indices = [[0, 1, 2, 3, 4, 5]]
    pharmacophore._aromatic_pharmacophoric_points(
        lig_centers, lig_indices,
        rec_centers, rec_indices,
        frame=0
    )
    assert len(pharmacophore[0]) == 0


def test_aromatic_pharmacophoric_points_pstack(mocker, ligand_chem_feats, receptor_chem_feats):
    pharmacophore = empty_pharmacophore_one_frame()
    pharmacophore._pl_complex = mocker.Mock()
    pharmacophore._pl_complex.coords = puw.quantity(np.array([
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
    ]), "angstroms")
    pharmacophore._pl_complex.coords = pharmacophore._pl_complex.coords.reshape((1, 12, 3))
    assert pharmacophore._pl_complex.coords.shape == (1, 12, 3)
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
    assert np.allclose(pharmacophore[0][0].direction,
                       np.array([1/np.sqrt(10), 0, 3/np.sqrt(10)]))


def test_aromatic_pharmacophoric_points_exceeds_max_offset(mocker, ligand_chem_feats, receptor_chem_feats):
    pharmacophore = empty_pharmacophore_one_frame()
    pharmacophore._pl_complex = mocker.Mock()
    pharmacophore._pl_complex.coords = puw.quantity(np.array([
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
    ]), "angstroms")
    pharmacophore._pl_complex.coords = pharmacophore._pl_complex.coords.reshape((1, 12, 3))
    assert pharmacophore._pl_complex.coords.shape == (1, 12, 3)
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
    pharmacophore._pl_complex.coords = puw.quantity(np.array([
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
    ]), "angstroms")
    pharmacophore._pl_complex.coords = pharmacophore._pl_complex.coords.reshape((1, 12, 3))
    assert pharmacophore._pl_complex.coords.shape == (1, 12, 3)
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
    assert np.allclose(pharmacophore[0][0].direction,
                       np.array([1., 0., 0.]))


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


def setup_extract(mocker, ligand_chem_feats, receptor_chem_feats, index=None):
    """ Set up a test for the extract method of LigandReceptorPharmacophore. """
    pharma = LigandReceptorPharmacophore()
    pharma._pl_complex = mocker.Mock()
    pl = pharma._pl_complex

    ligand_side_effect = [
        (ligand_chem_feats.hyd_cent, []),
        (ligand_chem_feats.charge_cent, []),
        (ligand_chem_feats.charge_cent, []),
        (ligand_chem_feats.aro_cent, ligand_chem_feats.aro_ind)
    ]
    receptor_side_effect = [
        (receptor_chem_feats.hyd_cent, []),
        (receptor_chem_feats.charge_cent, []),
        (receptor_chem_feats.charge_cent, []),
        (receptor_chem_feats.aro_cent, receptor_chem_feats.aro_ind)
    ]

    if index is None:
        pl.ligand_features.side_effect = ligand_side_effect
        pl.receptor_features.side_effect = receptor_side_effect
    else:
        # Take just one feature
        pl.ligand_features.side_effect = ligand_side_effect[index:index + 1]
        pl.receptor_features.side_effect = receptor_side_effect[index:index + 1]

    pl.coords = puw.quantity(np.array([
        # Aromatic ring ligand
        [1, 0, 0],
        [1 / 2, np.sqrt(3) / 2, 0],
        [-1 / 2, np.sqrt(3) / 2, 0],
        [-1, 0, 0],
        [-1 / 2, -np.sqrt(3) / 2, 0],
        [1 / 2, -np.sqrt(3) / 2, 0],
        # Aromatic ring receptor
        [2, 0, 3],
        [3 / 2, np.sqrt(3) / 2, 3],
        [1 / 2, np.sqrt(3) / 2, 3],
        [0, 0, 3],
        [1 / 2, -np.sqrt(3) / 2, 3],
        [3 / 2, -np.sqrt(3) / 2, 3],
        # Hydrogen bonding atoms receptor
        [0, 0, 0],  # donor atom
        [1, 1, 1],  # hydrogen atom
        [2, 2, 2],  # acceptor atom
        # Hydrogen bonding atoms ligand
        [3, 3, 3],  # acceptor atom
        [4, 4, 4],  # hydrogen atom
        [5, 5, 5],  # donor atom
    ]), "angstroms")
    # Coords must be a 3D array
    pl.coords = pl.coords.reshape((1, 18, 3))
    assert pl.coords.shape == (1, 18, 3)
    pl.hbond_indices.return_value = np.array([
        [12, 13, 14],
        [15, 16, 17],
    ])
    pl.lig_indices = [0, 1, 2, 3, 4, 5, 14, 15, 16]

    return pharma, pl


def test_extract_all_features(
        mocker, ligand_chem_feats, receptor_chem_feats):
    pharma, pl = setup_extract(mocker, ligand_chem_feats, receptor_chem_feats)
    pharma.extract("EST:B")

    pl.prepare.assert_called_once_with(
        lig_id="EST:B", smiles="", add_hydrogens=True
    )

    assert len(pharma[0]) == 6
    expected_features = {
        "hb acceptor",
        "aromatic ring",
        "hb donor",
        "hydrophobicity",
        "positive charge",
        "negative charge",
    }
    actual_features = set()
    for pnt in pharma[0]:
        actual_features.add(pnt.feature_name)
    assert actual_features == expected_features


def setup_test_extract_feature(mocker, lig_feats,
                               rec_feats, feature, index):
    """ Set up a test for the extract method with a particular feature."""
    pharma, pl = setup_extract(
        mocker, lig_feats, rec_feats, index)
    pharma.extract("EST:B", features=[feature])
    assert len(pharma[0]) == 1
    assert pharma[0][0].feature_name == feature


def test_extract_hydrophobic_features(
        mocker, ligand_chem_feats, receptor_chem_feats):
    setup_test_extract_feature(
        mocker, ligand_chem_feats,
        receptor_chem_feats, "hydrophobicity", 0
    )


def test_extract_positive_charge_features(
        mocker, ligand_chem_feats, receptor_chem_feats):
    setup_test_extract_feature(
        mocker, ligand_chem_feats,
        receptor_chem_feats, "positive charge", 1
    )


def test_extract_negative_charge_features(
        mocker, ligand_chem_feats, receptor_chem_feats):
    setup_test_extract_feature(
        mocker, ligand_chem_feats,
        receptor_chem_feats, "negative charge", 2
    )


def test_extract_aromatic_features(
        mocker, ligand_chem_feats, receptor_chem_feats):
    setup_test_extract_feature(
        mocker, ligand_chem_feats,
        receptor_chem_feats, "aromatic ring", 3
    )


def test_extract_acceptor_features(
        mocker, ligand_chem_feats, receptor_chem_feats):
    setup_test_extract_feature(
        mocker, ligand_chem_feats,
        receptor_chem_feats, "hb acceptor", None
    )


def test_extract_donor_features(
        mocker, ligand_chem_feats, receptor_chem_feats):
    setup_test_extract_feature(
        mocker, ligand_chem_feats,
        receptor_chem_feats, "hb donor", None
    )


def test_adding_single_frame_to_view_updates_components(pharmacophore_one_frame):
    view = nv.NGLWidget()
    assert len(view._ngl_component_ids) == 0
    pharmacophore_one_frame.add_to_view(view, 0)
    assert len(view._ngl_component_ids) == 3


@pytest.fixture()
def pharma_with_pl_complex(small_pdb_with_ligand):
    """ Returns a pharmacophore with a PLComplex object.
    """
    ph = LigandReceptorPharmacophore()
    ph.load_receptor(small_pdb_with_ligand)
    ph.add_frame()

    acceptor = PharmacophoricPoint(
        "hb acceptor",
        puw.quantity([1.0, 1.0, 1.0], "angstroms"),
        puw.quantity(1.0, "angstroms")

    )

    ph.add_points_to_frame([acceptor], 0)

    return ph


def test_show_all_single_frame(mocker, pharma_with_pl_complex):
    mock_nv = mocker.patch(
        lr_module + ".nv"
    )
    pharma = pharma_with_pl_complex
    pharma._pl_complex._receptor_indices = list(range(0, 10))

    view = pharma.show(
        frame=0, ligand=True, receptor=True, points=True
    )

    mock_nv.show_mdtraj.assert_called_once()
    mock_nv.show_rdkit.assert_not_called()
    _, args, _ = mock_nv.show_mdtraj.mock_calls[0]
    assert args[0].n_atoms == 10

    view.add_component.assert_called_once()


def test_show_ligand_no_receptor(mocker, pharmacophore_one_frame):
    mock_nv = mocker.patch(
        lr_module + ".nv"
    )
    pharma = pharmacophore_one_frame
    pharma._pl_complex = mocker.Mock()
    view = pharma.show(
        frame=0, ligand=True, receptor=False, points=True
    )

    mock_nv.show_rdkit.assert_called_once()
    mock_nv.show_mdtraj.assert_not_called()
    view.add_component.assert_not_called()


def test_show_receptor_no_ligand(mocker, pharma_with_pl_complex):
    mock_nv = mocker.patch(
        lr_module + ".nv"
    )
    pharma = pharma_with_pl_complex
    pharma._pl_complex._receptor_indices = list(range(0, 10))
    view = pharma.show(
        frame=0, ligand=False, receptor=True, points=True
    )

    mock_nv.show_mdtraj.assert_called_once()
    mock_nv.show_rdkit.assert_not_called()
    _, args, _ = mock_nv.show_mdtraj.mock_calls[0]
    assert args[0].n_atoms == 10

    view.add_component.assert_not_called()


def test_show_custom_indices_and_ligand(mocker, pharmacophore_one_frame):
    mocker.patch(
        lr_module + ".nv"
    )
    mock_add = mocker.patch(
        lr_class + ".add_to_view"
    )
    pharma = pharmacophore_one_frame
    pharma._pl_complex = mocker.Mock()

    view = pharma.show(
        frame=0, ligand=True, receptor=True, points=True,
        indices=np.array([0, 1, 2])
    )

    pharma.receptor.slice_traj.assert_called_once()
    _, args, _ = pharma.receptor.slice_traj.mock_calls[0]
    assert np.all(args[0] == np.array([0, 1, 2]))
    assert args[1] == 0  # Frame

    view.add_component.assert_called_once()

    mock_add.assert_called_once_with(view, 0)
    assert len(view.representations) == 1
    assert view.representations[0]["type"] == "ball+stick"
    assert view.representations[0]["params"]["sele"] == "all"


def test_show_all_traj_contains_multiple_frames(mocker, small_trajectory_path):
    mock_nv = mocker.patch(
        lr_module + ".nv"
    )

    ph = LigandReceptorPharmacophore()
    ph.load_receptor(small_trajectory_path)
    ph.receptor._receptor_indices = list(range(10))
    ph.receptor.get_lig_conformer = mocker.Mock()

    ph.add_frame()
    ph.add_frame()
    acceptor = PharmacophoricPoint(
        "hb acceptor",
        puw.quantity([1.0, 1.0, 1.0], "angstroms"),
        puw.quantity(1.0, "angstroms")

    )
    ph.add_points_to_frame([acceptor], 1)

    view = ph.show(
        frame=1, ligand=True, receptor=True, points=True
    )

    mock_nv.show_mdtraj.assert_called_once()
    mock_nv.show_rdkit.assert_not_called()
    _, args, _ = mock_nv.show_mdtraj.mock_calls[0]
    traj = args[0]
    assert traj.n_atoms == 10
    assert traj.n_frames == 1

    ph.receptor.get_lig_conformer.assert_called_once_with(1)
    view.add_component.assert_called_once()
