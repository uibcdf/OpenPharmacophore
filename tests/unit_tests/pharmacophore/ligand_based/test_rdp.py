import openpharmacophore.pharmacophore.ligand_based.rdp as rdp
import numpy as np
import pytest


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

    assert len(boxes) == 1
    assert len(boxes[0]) == 5
    assert boxes[0][0].id == (0, 2)
    assert boxes[0][1].id == (0, 3)
    assert boxes[0][2].id == (1, 0)
    assert boxes[0][3].id == (2, 2)
    assert boxes[0][4].id == (3, 0)


def test_retrieve_cps():
    assert False
