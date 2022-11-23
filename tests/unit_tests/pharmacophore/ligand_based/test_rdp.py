import openpharmacophore.pharmacophore.ligand_based.rdp as rdp
import numpy as np
import pytest
from unittest.mock import Mock


def test_init_feature_list():
    f_list = rdp.FeatureList("AAR", (0, 1), np.array([1.0, 1.0, 1.0]))
    assert f_list.id == (0, 1)
    assert f_list.variant == "AAR"
    assert np.all(f_list.distances == np.array([1.0, 1.0, 1.0]))


def test_init_feat_list_distances_shape_incorrect_raises_error():
    with pytest.raises(ValueError):
        rdp.FeatureList("AAR", (0, 1), np.array([1.0, 1.0]))


def test_init_flcontainer():
    rdp.FLContainer()  # Should not raise
    rdp.FLContainer(bin=(0, 1))
    rdp.FLContainer(variant="")


def test_flcontainer_append():
    container = rdp.FLContainer()
    container.append(rdp.FeatureList("AAR", (0, 0), np.array([4.8, 2.8, 7.0])))
    assert container.variant == "AAR"
    assert container.mols == {0}
    assert len(container.flists) == 1
    assert container[0].id == (0, 0)

    container.append(rdp.FeatureList("AAR", (1, 0), np.array([4.8, 2.8, 7.0])))
    assert container.mols == {0, 1}
    assert len(container.flists) == 2
    assert container[1].id == (1, 0)


def test_flcontainer_append_different_variant_raises_error():
    container = rdp.FLContainer()
    container.append(rdp.FeatureList("AAR", (0, 0), np.array([4.8, 2.8, 7.0])))
    with pytest.raises(ValueError):
        container.append(rdp.FeatureList("DHR", (0, 0), np.array([4.8, 2.8, 7.0])))


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


def test_recursive_partitioning():

    lists = [
        rdp.FeatureList("AAR", (0, 0), np.array([4.8, 2.8, 7.0])),
        rdp.FeatureList("AAR", (0, 1), np.array([3.1, 5.9, 3.0])),
        rdp.FeatureList("AAR", (0, 2), np.array([4.3, 3.7, 6.8])),
        rdp.FeatureList("AAR", (0, 3), np.array([3.7, 4.8, 6.7])),
        rdp.FeatureList("AAR", (0, 4), np.array([2.6, 4.2, 6.4])),
        rdp.FeatureList("AAR", (1, 0), np.array([4.4, 4.7, 6.1])),
        rdp.FeatureList("AAR", (1, 1), np.array([2.6, 5.5, 5.1])),
        rdp.FeatureList("AAR", (1, 2), np.array([5.5, 2.1, 4.9])),
        rdp.FeatureList("AAR", (1, 3), np.array([6.0, 5.7, 4.3])),
        rdp.FeatureList("AAR", (2, 0), np.array([4.9, 2.7, 2.9])),
        rdp.FeatureList("AAR", (2, 1), np.array([6.8, 5.6, 5.2])),
        rdp.FeatureList("AAR", (2, 2), np.array([4.5, 4.2, 5.8])),
        rdp.FeatureList("AAR", (3, 0), np.array([3.8, 5.2, 5.9])),
        rdp.FeatureList("AAR", (3, 1), np.array([5.1, 5.4, 4.6])),
    ]

    f_lists = rdp.FLContainer()
    for fl in lists:
        f_lists.append(fl)

    n_pair = f_lists[0].n_pairs
    assert n_pair == 3

    boxes = []
    rdp.recursive_partitioning(f_lists, 0, n_pair, boxes, 4)

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


def test_score_common_pharmacophores():
    box = rdp.FLContainer(variant="AAR")
    box.append(
        rdp.FeatureList("AAR", (0, 0), np.array([1.4, 0.6, 2.7]))
    )
    box.append(
        rdp.FeatureList("AAR", (1, 0), np.array([1.2, 0.8, 3.1]))
    )
    box.append(
        rdp.FeatureList("AAR", (2, 0), np.array([1.1, 0.9, 3.0]))
    )

    scores = rdp.score_common_pharmacophores(box)
    assert len(scores) == 3
    assert scores[0] == (0.9166666666666667, 1, 2)
    assert scores[1] == (0.7642977396044842, 0, 1)
    assert scores[2] == (0.75, 0, 2)


def test_score_common_pharmacophores_rmsd_cutoff_exceeded():
    box = rdp.FLContainer(variant="AAR")
    box.append(
        rdp.FeatureList("AAR", (0, 0), np.array([1.4, 0.6, 2.7]))
    )
    box.append(
        rdp.FeatureList("AAR", (1, 0), np.array([1.8, 1.0, 3.2]))
    )
    box.append(
        rdp.FeatureList("AAR", (2, 0), np.array([3.0, 1.6, 3.7],))
    )
    scores = rdp.score_common_pharmacophores(box)
    assert len(scores) == 1
    assert scores[0] == (0.311133512909042, 1, 2)


@pytest.fixture()
def ligands():
    """ Returns a list of rdp.Ligands with variants assigned.
    """
    mock_mol = Mock()
    ligands = [rdp.Ligands(mock_mol)] * 4
    for lig in ligands:
        lig.n_confs = 2

    variants = ["AAHP", "AAPR", "AADP", "AADPR"]
    for ii in range(len(ligands)):
        ligands[ii].variant = variants[ii]
    return ligands


def test_common_k_point_variants():
    variants = ["AAHP", "AAPR", "AADP", "AADPR"]
    common_variants = rdp.common_k_point_variants(
        variants, n_points=3, min_actives=4)
    assert common_variants == ["AAP"]


def test_common_k_point_variants_min_actives_less_than_variants():
    variants = ["AAHP", "AAPR", "AADP", "AADPR"]
    common_variants = rdp.common_k_point_variants(
        variants, n_points=3, min_actives=2)
    assert common_variants == ["AAP", "AAR", "APR", "AAD", "ADP"]


def test_common_k_point_feature_lists(mocker, ligands):
    mocker.patch("openpharmacophore.pharmacophore.ligand_based.rdp.Ligand.distances",
                 return_value=np.array([1., 1., 1.]))
    mock_mol = mocker.Mock()
    mock_mol.GetNumConformers.return_value = 2

    k_variants = ["AAP", "AAR", "APR", "AAD", "ADP"]
    containers = rdp.common_k_point_feature_lists(ligands, k_variants)
    assert len(containers) == 5

    assert containers[0].variant == "AAP"
    assert len(containers[0]) == 8

    assert containers[1].variant == "AAR"
    assert len(containers[1]) == 4

    assert containers[2].variant == "APR"
    assert len(containers[2]) == 4

    assert containers[3].variant == "AAD"
    assert len(containers[3]) == 4

    assert containers[4].variant == "ADP"
    assert len(containers[4]) == 4


def test_find_common_pharmacophores():
    pass
