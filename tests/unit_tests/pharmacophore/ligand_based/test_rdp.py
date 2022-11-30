import openpharmacophore.pharmacophore.ligand_based.rdp as rdp
import numpy as np
import pytest
import pyunitwizard as puw
from unittest.mock import Mock
from collections import Counter


# Feature List Class Tests


def test_init_feature_list():
    f_list = rdp.FeatureList("AAR", (0, 1, 2), (0, 1),
                             np.array([1.0, 1.0, 1.0]))
    assert f_list.id == (0, 1)
    assert f_list.variant == "AAR"
    assert f_list.var_ind == (0, 1, 2)
    assert np.all(f_list.distances == np.array([1.0, 1.0, 1.0]))


def test_init_feat_list_distances_shape_incorrect_raises_error():
    with pytest.raises(ValueError):
        rdp.FeatureList("AAR", (0, 1, 2), (0, 1), np.array([1.0, 1.0]))


# FLContainer Class Tests

def test_init_flcontainer():
    rdp.FLContainer()  # Should not raise
    rdp.FLContainer(bin=(0, 1))
    rdp.FLContainer(variant="")


def test_flcontainer_append():
    container = rdp.FLContainer()
    container.append(rdp.FeatureList("AAR", (0, 1, 2), (0, 0),
                                     np.array([4.8, 2.8, 7.0])))
    assert container.variant == "AAR"
    assert container.mols == {0}
    assert len(container.flists) == 1
    assert container[0].id == (0, 0)

    container.append(rdp.FeatureList("AAR", (0, 1, 2),
                                     (1, 0), np.array([4.8, 2.8, 7.0])))
    assert container.mols == {0, 1}
    assert len(container.flists) == 2
    assert container[1].id == (1, 0)


def test_flcontainer_append_different_variant_raises_error():
    container = rdp.FLContainer()
    container.append(rdp.FeatureList("AAR", (0, 1, 2),
                                     (0, 0), np.array([4.8, 2.8, 7.0])))
    with pytest.raises(ValueError):
        container.append(rdp.FeatureList("DHR", (0, 1, 2),
                                         (0, 0), np.array([4.8, 2.8, 7.0])))


# Ligand Class tests

def test_init_ligand(mocker):
    lig = rdp.Ligand(mocker.Mock, {"A": [(1,)]})  # Should not raise


def test_update_variant(mocker):
    feats = {
        "R": [(1, 2, 3)],
        "H": [(4,)],
        "N": [(5,)],
        "P": [],
        "A": [(6,), (7,)],
        "D": [],
    }
    mock_mol = mocker.Mock()
    mock_mol.GetNumConformers.return_value = 2
    lig = rdp.Ligand(mock_mol, feats)
    lig._update_variant()
    assert lig.variant == "AAHNR"
    assert lig.feat_count == {
        "A": 2,
        "H": 1,
        "N": 1,
        "R": 1,
    }


def test_k_distances_values_precomputed():
    lig = rdp.Ligand(None, {})
    lig.variant = "AAHR"
    lig.distances = np.array([
        [[-1., 1., 2., 3.],
         [1., -1., 4., 5.],
         [2., 4., -1., 6.],
         [3., 5., 6., -1.],
         ],
    ])
    assert lig.distances.shape == (1, 4, 4)

    expected = np.array([1., 3., 5.])
    assert np.all(lig.k_distances((0, 1, 3), 0) == expected)

    expected = np.array([4., 5., 6.])
    assert np.all(lig.k_distances((1, 2, 3), 0) == expected)


def test_interpoint_distances(mocker):
    mocker.patch(
        "openpharmacophore.pharmacophore.ligand_based.rdp.feature_centroids",
        side_effect=[
            np.array([0., 0., 0.]),
            np.array([1., 1., 1.]),
            np.array([3., 3., 3.]),
        ]
    )
    feats = {
        "A": [(0,), (6,)],
        "R": [(1, 2, 3)],
    }
    lig = rdp.Ligand(None, feats)
    lig.variant = "AAR"
    lig.feat_count = Counter(lig.variant)
    lig.distances = np.array([
        [[-1, -1, -1],
         [-1, -1, -1],
         [-1, -1, -1]],
    ], dtype=float)
    lig.interpoint_distances(0)

    d1 = np.sqrt(3)
    d2 = np.sqrt(27)
    d3 = 2 * d1
    expected = np.array([
        [[-1, d1, d2],
         [d1, -1, d3],
         [d2, d3, -1]],
    ])
    assert np.allclose(lig.distances, expected)


# Recursive partitioning and common pharmacophores tests


def test_nearest_bins():
    bin_size = 1.0
    bins = np.arange(0, 16, step=bin_size)

    nearest = rdp.nearest_bins(2.3, bin_size)
    assert nearest == (1, 2)
    assert bins[nearest[0]] == 1
    assert bins[nearest[1]] == 2

    nearest = rdp.nearest_bins(4.9, bin_size)
    assert nearest == (4, 5)
    assert bins[nearest[0]] == 4
    assert bins[nearest[1]] == 5

    bin_size = 2.0
    bins = np.arange(0, 16, step=bin_size)

    nearest = rdp.nearest_bins(2.5, bin_size)
    assert nearest == (0, 1)
    assert bins[nearest[0]] == 0
    assert bins[nearest[1]] == 2

    nearest = rdp.nearest_bins(3.9, bin_size)
    assert nearest == (2, 3)
    assert bins[nearest[0]] == 4
    assert bins[nearest[1]] == 6


@pytest.fixture()
def fl_container():
    lists = [
        rdp.FeatureList("AAR", (0, 1, 2), (0, 0), np.array([4.8, 2.8, 7.0])),
        rdp.FeatureList("AAR", (0, 1, 2), (0, 1), np.array([3.1, 5.9, 3.0])),
        rdp.FeatureList("AAR", (0, 1, 2), (0, 2), np.array([4.3, 3.7, 6.8])),
        rdp.FeatureList("AAR", (0, 1, 2), (0, 3), np.array([3.7, 4.8, 6.7])),
        rdp.FeatureList("AAR", (0, 1, 2), (0, 4), np.array([2.6, 4.2, 6.4])),
        rdp.FeatureList("AAR", (0, 1, 2), (1, 0), np.array([4.4, 4.7, 6.1])),
        rdp.FeatureList("AAR", (0, 1, 2), (1, 1), np.array([2.6, 5.5, 5.1])),
        rdp.FeatureList("AAR", (0, 1, 2), (1, 2), np.array([5.5, 2.1, 4.9])),
        rdp.FeatureList("AAR", (0, 1, 2), (1, 3), np.array([6.0, 5.7, 4.3])),
        rdp.FeatureList("AAR", (0, 1, 2), (2, 0), np.array([4.9, 2.7, 2.9])),
        rdp.FeatureList("AAR", (0, 1, 2), (2, 1), np.array([6.8, 5.6, 5.2])),
        rdp.FeatureList("AAR", (0, 1, 2), (2, 2), np.array([4.5, 4.2, 5.8])),
        rdp.FeatureList("AAR", (0, 1, 2), (3, 0), np.array([3.8, 5.2, 5.9])),
        rdp.FeatureList("AAR", (0, 1, 2), (3, 1), np.array([5.1, 5.4, 4.6])),
    ]
    f_lists = rdp.FLContainer()
    for fl in lists:
        f_lists.append(fl)
    return f_lists


def test_pharmacophore_partitioning(fl_container):
    boxes = rdp.pharmacophore_partitioning(fl_container, 4)
    assert len(boxes) == 3

    ids = [b.id for b in boxes[0]]
    assert ids == [
        (0, 4), (1, 0), (1, 1), (2, 2), (3, 0)
    ]
    ids = [b.id for b in boxes[1]]
    assert ids == [
        (0, 2), (0, 3), (0, 4), (1, 0), (2, 2), (3, 0)
    ]
    ids = [b.id for b in boxes[2]]
    assert ids == [
        (0, 2), (0, 3), (1, 0), (2, 2), (3, 0)
    ]


def test_pharmacophore_partitioning_min_actives_less_than_ligands_size(fl_container):
    boxes = rdp.pharmacophore_partitioning(fl_container, 3)
    assert len(boxes) == 6

    ids = [b.id for b in boxes[2]]
    assert ids == [(0, 3), (1, 0), (3, 0)]
    ids = [b.id for b in boxes[3]]
    assert ids == [(1, 0), (2, 2), (3, 0), (3, 1)]
    ids = [b.id for b in boxes[5]]
    assert ids == [(0, 3), (1, 0), (3, 0)]


def test_score_common_pharmacophores():
    box = rdp.FLContainer(variant="AAR")
    box.append(
        rdp.FeatureList("AAR", (0, 1, 2), (0, 0), np.array([1.4, 0.6, 2.7]))
    )
    box.append(
        rdp.FeatureList("AAR", (0, 1, 2), (1, 0), np.array([1.2, 0.8, 3.1]))
    )
    box.append(
        rdp.FeatureList("AAR", (0, 1, 2), (2, 0), np.array([1.1, 0.9, 3.0]))
    )

    scores = rdp.score_common_pharmacophores(box)
    assert len(scores) == 3
    assert scores[0] == (0.9166666666666667, 1, 2)
    assert scores[1] == (0.7642977396044842, 0, 1)
    assert scores[2] == (0.75, 0, 2)


def test_score_common_pharmacophores_rmsd_cutoff_exceeded():
    box = rdp.FLContainer(variant="AAR")
    box.append(
        rdp.FeatureList("AAR", (0, 1, 2), (0, 0), np.array([1.4, 0.6, 2.7]))
    )
    box.append(
        rdp.FeatureList("AAR", (0, 1, 2), (1, 0), np.array([1.8, 1.0, 3.2]))
    )
    box.append(
        rdp.FeatureList("AAR", (0, 1, 2), (2, 0), np.array([3.0, 1.6, 3.7], ))
    )
    scores = rdp.score_common_pharmacophores(box)
    assert len(scores) == 1
    assert scores[0] == (0.311133512909042, 1, 2)


@pytest.fixture()
def ligands():
    """ Returns a list of rdp.Ligands with variants assigned.
    """
    mock_mol = Mock()
    mock_mol.GetNumConformers.return_value = 2
    ligands = [rdp.Ligand(mock_mol, {}) for _ in range(4)]

    variants = ["AAHP", "AAPR", "AADP", "AADPR"]
    for ii in range(len(ligands)):
        ligands[ii].variant = variants[ii]
        ligands[ii].feat_count = Counter(variants[ii])
    return ligands


def test_common_k_point_feat_lists(ligands):
    common_variants = rdp.common_k_point_variants(
        ligands, n_points=3, min_actives=4)

    assert len(common_variants) == 4

    assert common_variants[0].name == "AAP"
    assert common_variants[0].indices == (0, 1, 3)
    assert common_variants[0].mol == 0

    assert common_variants[1].name == "AAP"
    assert common_variants[1].indices == (0, 1, 2)
    assert common_variants[1].mol == 1

    assert common_variants[2].name == "AAP"
    assert common_variants[2].indices == (0, 1, 3)
    assert common_variants[2].mol == 2

    assert common_variants[3].name == "AAP"
    assert common_variants[3].indices == (0, 1, 3)
    assert common_variants[3].mol == 3


def test_common_k_point_feat_lists_min_actives_less_than_variants(ligands):
    common = rdp.common_k_point_variants(
        ligands, n_points=3, min_actives=2)

    assert len(common) == 16
    expected_var = ["AAP"] * 4
    expected_var += ["AAR"] * 2
    expected_var += ["APR"] * 4
    expected_var += ["AAD"] * 2
    expected_var += ["ADP"] * 4
    var = [f.name for f in common]
    assert var == expected_var

    expected_mols = [0, 1, 2, 3, 1, 3, 1, 1, 3, 3, 2, 3, 2, 2, 3, 3]
    mols = [f.mol for f in common]
    assert mols == expected_mols


def test_common_k_point_feature_lists(mocker, ligands):
    mocker.patch("openpharmacophore.pharmacophore.ligand_based.rdp.Ligand.k_distances",
                 return_value=np.array([1., 1., 1.]))

    k_variants = rdp.common_k_point_variants(ligands, n_points=3, min_actives=2)
    containers = rdp.common_k_point_feature_lists(ligands, k_variants)

    assert len(containers) == 5
    expected_variants = ["AAP", "AAR", "APR", "AAD", "ADP"]
    expected_size = [8, 4, 8, 4, 8]
    expected_mols = [{0, 1, 2, 3}, {1, 3}, {1, 3},
                     {2, 3}, {2, 3}]
    var_names = [c.variant for c in containers]
    size = [len(c) for c in containers]
    mols = [c.mols for c in containers]

    assert var_names == expected_variants
    assert size == expected_size
    assert mols == expected_mols


def test_feat_list_to_pharma(mocker):
    mock_centroids = mocker.patch(
        "openpharmacophore.pharmacophore.ligand_based.rdp.feature_centroids",
        side_effect=[
            np.array([0., 0., 0.]),
            np.array([1., 1., 1.]),
            np.array([2., 2., 2.]),
        ]
    )
    feats = {
        "A": [(9,), (10,)],
        "P": [(2, 3, 4)],
        "R": [(6, 7, 8)],
    }
    ligand = rdp.Ligand(None, feats)
    ligand.variant = "AAAPR"

    fl = rdp.FeatureList("APR", (1, 3, 4), (0, 2), np.zeros((3,)))
    pharma = rdp.feat_list_to_pharma(fl, ligand)

    assert mock_centroids.call_count == 3
    calls = [
        mocker.call(None, 2, (10, )),
        mocker.call(None, 2, (2, 3, 4)),
        mocker.call(None, 2, (6, 7, 8)),
    ]
    assert mock_centroids.call_args_list == calls

    assert len(pharma) == 3
    feat_names = [p.feature_name for p in pharma]
    assert feat_names == ["hb acceptor", "positive charge", "aromatic ring"]
    assert np.all(puw.get_value(pharma[0].center) == np.zeros((3,)))
    assert np.all(puw.get_value(pharma[1].center) == np.ones((3,)))
    assert np.all(puw.get_value(pharma[2].center) == np.ones((3,)) * 2)
