import openpharmacophore.pharmacophore.pl_interactions as pli
import openpharmacophore.data as data
import numpy as np
import pytest
import mdtraj as mdt


@pytest.fixture
def trajectory():
    return mdt.load(data.pdb["1qku"])


def test_find_ligands(trajectory):
    ligands = pli.find_ligands_in_traj(trajectory)
    assert len(ligands) == 3
    assert ligands == [
        "EST:D",
        "EST:E",
        "EST:F",
    ]


def test_get_ligand_atom_indices(trajectory):
    traj = trajectory
    ligand_indices = pli.get_ligand_atom_indices(traj, "EST:D")
    expected_indices = list(range(5940, 5960))
    assert ligand_indices == expected_indices


def test_ligand_centroid(trajectory):
    traj = trajectory
    ligand_id = "EST:D"
    centroid = pli.ligand_centroid(traj, ligand_id)

    coordinates = np.array(
        [[[10.4106, 1.7203, 2.4775],
          [10.2995, 1.7834, 2.537],
          [10.1695, 1.7355, 2.512],
          [10.0598, 1.799, 2.5704],
          [10.1506, 1.624, 2.4274],
          [10.2621, 1.5588, 2.366],
          [10.2371, 1.4379, 2.2735],
          [10.3644, 1.3753, 2.2086],
          [10.4898, 1.3873, 2.2953],
          [10.5178, 1.5388, 2.3261],
          [10.3957, 1.6078, 2.3918],
          [10.6462, 1.5459, 2.4125],
          [10.7711, 1.4803, 2.3508],
          [10.7463, 1.3343, 2.3124],
          [10.617, 1.327, 2.2242],
          [10.6228, 1.1821, 2.1792],
          [10.7701, 1.1713, 2.1263],
          [10.8494, 1.2719, 2.2135],
          [10.961, 1.2027, 2.2746],
          [10.7379, 1.2449, 2.4419]]],
        dtype=np.float32)
    expected_centroid = np.mean(coordinates, axis=1)[0]
    assert centroid.shape == (3, )
    assert np.allclose(centroid, expected_centroid)


def test_maximum_distance():
    centroid = np.array([0., 0., 0.])
    coordinates = np.array([
        [1., 1., 1.],
        [2., 2., 2.],
        [3., 3., 3.],
        [4., 4., 4.],
    ])
    assert pli.maximum_distance(centroid, coordinates) == np.sqrt(48)


def test_maximum_ligand_extent(trajectory):
    assert np.allclose(pli.ligand_maximum_extent(trajectory, "EST:D"), 0.59849)


def test_get_binding_site_atoms_indices():
    ligand_maximum_extent = 0.15
    ligand_centroid = np.array([0.0, 0.0, 0.0])
    coordinates = np.array([
        [0.01, 0.01, 0.01],
        [0.02, 0.02, 0.02],
        [0.3, 0.3, 0.3],
        [0.4, 0.4, 0.4],
        [0.5, 0.5, 0.5],
        [1., 1., 1.],
        [2., 2., 2.],
    ])
    bsite_indices = pli.get_binding_site_atoms_indices(
        ligand_centroid, ligand_maximum_extent, coordinates)
    assert isinstance(bsite_indices, np.ndarray)

    expected_indices = np.array([2, 3, 4])
    assert np.all(bsite_indices == expected_indices)


@pytest.fixture()
def smart_patterns_mock():
    return {
        'N=[CX3](N)-N': 'PosIonizable',
        '[#16!H0]': 'Donor',
        '[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]': 'Acceptor',
        '[C&r3]1~[C&r3]~[C&r3]1': 'Hydrophobe',
        'a1aaaa1': 'Aromatic',
        'c1nn[nH1]n1': 'NegIonizable'}


side_effect = [
        ((0, 1), (1, 2)),
        ((2, 3), ),
        ((3, 4), (5, 6)),
        ((7, 8, 9),),
        ((10, 11, 12), (13, 14, 15)),
        (),
    ]


def test_ligand_chemical_features(mocker, smart_patterns_mock):
    mocker.patch.dict(
        pli.smarts_patterns,
        smart_patterns_mock,
        clear=True
    )
    mock_mol = mocker.Mock()
    mock_mol.GetSubstructMatches.side_effect = side_effect
    features = pli.chemical_features(mock_mol, None)
    expected_features = {
        'PosIonizable': [(0, 1), (1, 2)],
        'Donor': [(2, 3)],
        'Acceptor': [(3, 4), (5, 6)],
        'Hydrophobe': [(7, 8, 9)],
        'Aromatic': [(10, 11, 12), (13, 14, 15)]
    }
    assert features == expected_features


def test_binding_site_chemical_features(mocker, smart_patterns_mock):
    mocker.patch.dict(
        pli.smarts_patterns,
        smart_patterns_mock,
        clear=True
    )
    mock_mol = mocker.Mock()
    mock_mol.GetSubstructMatches.side_effect = side_effect
    bs_indices = np.array([5, 6, 7, 8, 9, 10, 11, 12])
    features = pli.chemical_features(mock_mol, bs_indices)
    expected_features = {
        "Acceptor": [(5, 6)],
        'Hydrophobe': [(7, 8, 9)],
        'Aromatic': [(10, 11, 12)]
    }
    assert features == expected_features


def test_features_centroid():
    coords = np.array([[
        [1., 1., 1.],
        [2., 2., 2.],
        [3., 3., 3.],
        [4., 4., 4.],
        [5., 5., 5.],
        [6., 6., 6.],
    ]])
    assert coords.shape == (1, 6, 3)
    features = {
        "Acceptor": [(1,)],
        "Hydrophobe": [(2, 3), (4,)]
    }
    expected_feats_centroid = {
        "Acceptor": [np.array([2., 2., 2.])],
        "Hydrophobe": [np.array([3.5, 3.5, 3.5]),
                       np.array([5., 5., 5.])]
    }
    feat_centroid = pli.features_centroid(features, coords)
    assert len(feat_centroid) == 2
    assert feat_centroid["Acceptor"][0].shape == (3,)
    assert feat_centroid["Hydrophobe"][0].shape == (3,)
    assert np.all(feat_centroid["Acceptor"][0] == expected_feats_centroid["Acceptor"][0])
    assert np.all(feat_centroid["Hydrophobe"][0] == expected_feats_centroid["Hydrophobe"][0])
    assert np.all(feat_centroid["Hydrophobe"][1] == expected_feats_centroid["Hydrophobe"][1])
