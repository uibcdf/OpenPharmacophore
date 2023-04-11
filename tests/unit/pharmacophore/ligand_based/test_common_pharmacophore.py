from collections import defaultdict
import pyunitwizard as puw
import numpy as np

from openpharmacophore import Pharmacophore, PharmacophoricPoint
from openpharmacophore.pharmacophore.align import align_pharmacophores
import openpharmacophore.pharmacophore.ligand_based.common_pharmacophore as cp
from openpharmacophore.molecular_systems.chem_feats import ChemFeat, ChemFeatContainer


class TestCommonPharmacophoreFinder:

    def test_init_common_pharmacophore_finder(self):
        scoring_fn_params = {
            "point_weight": 1.0,
            "rmsd_cutoff": puw.quantity(1.2, "angstroms"),
            "vector_weight": 1.0,
            "cos_cutoff": 0.5
        }
        min_dist = puw.quantity(1, "angstroms")
        max_dist = puw.quantity(10, "angstroms")
        bin_size = min_dist

        finder = cp.CommonPharmacophoreFinder(
            n_points=3,
            min_actives=4,
            max_pharmacophores=10,
            scoring_fn_params=scoring_fn_params,
            max_dist=max_dist,
            min_dist=min_dist,
            bin_size=bin_size,
        )

        assert finder.min_dist == 1
        assert finder.max_dist == 10
        assert finder.bin_size == 1

        assert finder.scoring_fn.point_weight == 1.0
        assert finder.scoring_fn.vector_weight == 1.0
        assert finder.scoring_fn.rmsd_cutoff == puw.quantity(1.2, "angstroms")
        assert finder.scoring_fn.cos_cutoff == 0.5

    def test_find_common_pharmacophores(self):
        coords_1 = puw.quantity(np.array([
            [4., 3., 0.],
            [4., 0., 0.],
            [0., 0., 0.],
            [-2., -1, 0.],
        ]), "angstroms")

        # Conformer 1 belongs to the first molecule. Conformer 2 and 3 to the second
        conf_1 = ChemFeatContainer([
            ChemFeat(coords_1[0], "hb acceptor"),
            ChemFeat(coords_1[1], "hb donor"),
            ChemFeat(coords_1[2], "hydrophobicity"),
            ChemFeat(coords_1[3], "positive charge"),
        ])

        coords = puw.quantity(np.array([
            [3.125, 3.903123749, 0.],
            [4., 0., 0.],
            [0., 0., 0.],
            [-2., -1., 0.],
        ]), "angstroms")

        # Conformer 1 belongs to the first molecule. Conformer 2 and 3 to the second
        conf_2 = ChemFeatContainer([
            ChemFeat(coords[0], "hb acceptor"),
            ChemFeat(coords[1], "hb donor"),
            ChemFeat(coords[2], "hydrophobicity"),
            ChemFeat(coords[3], "negative charge"),
        ])

        coords = puw.quantity(np.array([
            [6.3333333333, 2.98142397, 0.],
            [6., 0., 0.],
            [0., 0., 0.],
            [-2., -1, 0.],
        ]), "angstroms")

        # Conformer 1 belongs to the first molecule. Conformer 2 and 3 to the second
        conf_3 = ChemFeatContainer([
            ChemFeat(coords[0], "hb acceptor"),
            ChemFeat(coords[1], "hb donor"),
            ChemFeat(coords[2], "hydrophobicity"),
            ChemFeat(coords[3], "aromatic ring"),
        ])

        molecules = [[conf_1], [conf_2, conf_3]]
        cp_finder = cp.CommonPharmacophoreFinder(align="coords")
        common_pharmacophores = cp_finder(molecules, n_points=3)

        radius = puw.quantity(1.0, "angstroms")
        assert common_pharmacophores == [
            Pharmacophore([
                PharmacophoricPoint("hb acceptor", coords_1[0], radius),
                PharmacophoricPoint("hb donor", coords_1[1], radius),
                PharmacophoricPoint("hydrophobicity", coords_1[2], radius),
            ],
                ref_struct=0,
                ref_mol=0,
                score=0.6506998281828833
            ),
        ]

    def test_get_feat_lists(self):
        coords_1 = puw.quantity(np.array([
            [4., 0., 0.],
            [0., 0., 0.],
        ]), "angstroms")

        conformer_1 = ChemFeatContainer([
            ChemFeat(coords_1[0], "hb donor"),
            ChemFeat(coords_1[1], "aromatic ring"),
        ])

        coords_2 = puw.quantity(np.array([
            [6., 0., 0.],
            [0., 0., 0.],
        ]), "angstroms")

        conformer_2 = ChemFeatContainer([
            ChemFeat(coords_2[0], "hb donor"),
            ChemFeat(coords_2[1], "aromatic ring"),
        ])

        feat_lists = cp.CommonPharmacophoreFinder._get_feat_lists([
            [conformer_1, conformer_2]
        ])
        assert len(feat_lists) == 1
        assert len(feat_lists[0]) == 2

        expected_distances = np.array([
            [4., ],
            [6., ],
        ])
        expected_coords = np.stack((coords_1, coords_2))

        f_list = feat_lists[0]
        assert f_list.variant == "DR"
        assert np.all(f_list.distances == expected_distances)
        assert np.all(f_list.coords == expected_coords)

    def test_common_k_point_variants(self):
        feat_lists = [
            cp.FeatureList("AADH", None, None),
            cp.FeatureList("ADHP", None, None),
            cp.FeatureList("ADHN", None, None),
        ]
        common = cp.CommonPharmacophoreFinder._common_k_point_variants(feat_lists, 3, 3)
        expected = {
            "ADH": [
                cp.KVariant(var_ind=(0, 2, 3), mol=0),
                cp.KVariant(var_ind=(1, 2, 3), mol=0),
                cp.KVariant(var_ind=(0, 1, 2), mol=1),
                cp.KVariant(var_ind=(0, 1, 2), mol=2),
            ]
        }
        assert expected == common

    def test_variant_sublists(self):
        common_variants = {
            "ADP": [cp.KVariant((0, 1, 2), 0), cp.KVariant((0, 1, 3), 1)]
        }
        coords = np.ones((2, 4, 3))
        distances = np.array([
            [1, 2, 3, 4, 5, 6],
            [7, 8, 9, 10, 11, 12],
        ])
        feat_lists = [
            cp.FeatureList("ADPR", distances, coords),
            cp.FeatureList("ADHP", distances * 2, coords),
        ]
        sublists = cp.CommonPharmacophoreFinder._variant_sublists(feat_lists, common_variants)

        expected = defaultdict(list)
        expected["ADP"] = [
            cp.KSubList(np.array([1, 2, 4]), mol_id=(0, 0), feat_ind=[0, 1, 2]),
            cp.KSubList(np.array([7, 8, 10]), mol_id=(0, 1), feat_ind=[0, 1, 2]),
            cp.KSubList(np.array([2, 6, 10]), mol_id=(1, 0), feat_ind=[0, 1, 3]),
            cp.KSubList(np.array([14, 18, 22]), mol_id=(1, 1), feat_ind=[0, 1, 3]),
        ]

        assert expected == sublists

    def test_recursive_partitioning(self):
        sublists = [
            cp.KSubList(np.array([3, 5, 4]), mol_id=(0, 0), feat_ind=[]),
            cp.KSubList(np.array([4, 5, 4]), mol_id=(1, 0), feat_ind=[]),
            cp.KSubList(np.array([3, 7, 6]), mol_id=(1, 1), feat_ind=[]),
        ]
        cp_finder = cp.CommonPharmacophoreFinder()
        surviving_boxes = cp_finder._recursive_partitioning(sublists, min_actives=2)
        assert surviving_boxes == [
            [cp.KSubList(np.array([3, 5, 4]), mol_id=(0, 0), feat_ind=[]),
             cp.KSubList(np.array([4, 5, 4]), mol_id=(1, 0), feat_ind=[]),
             ],
        ]

    def test_box_top_representative(self):
        cp.KSubList.ID = 1
        cp_finder = cp.CommonPharmacophoreFinder(align="coords")
        surviving_box = [
            cp.KSubList(np.array([3, 5, 4]), mol_id=(0, 0), feat_ind=[0, 1, 2]),
            cp.KSubList(np.array([4, 5, 4]), mol_id=(1, 0), feat_ind=[0, 1, 2]),
            cp.KSubList(np.array([6, 5, 4]), mol_id=(1, 1), feat_ind=[0, 1, 2]),
        ]
        scores = {}
        coords_1 = puw.quantity(np.array([
            [[4., 3., 0.],
             [4., 0., 0.],
             [0., 0., 0.],
             ]]), "angstroms")
        coords_2 = puw.quantity(np.array([
            [[3.125, 3.903123749, 0.],
             [4., 0., 0.],
             [0., 0., 0.],
             ],
            [[6.3333333333, 2.98142397, 0.],
             [6., 0., 0.],
             [0., 0., 0.],
             ],
        ]), "angstroms")

        feature_lists = [
            cp.FeatureList("AAR", distances=np.array([3, 5, 4]), coords=coords_1),
            cp.FeatureList("AAR", distances=np.array([3, 5, 4]), coords=coords_2),
        ]

        top = cp_finder._box_top_representative(surviving_box, scores, feature_lists)
        assert top == surviving_box[0]
        assert top.score == 0.6506998281828833
        assert scores == {
            (1, 2): 0.6506998281828833,
            (1, 3): 0.14632288133597426,
        }

    def test_get_pharmacophores(self):
        queue = cp.FeatListQueue(size=2)
        sublist1 = cp.KSubList(np.ones(2), (0, 0), [0, 1], score=0.5)
        sublist2 = cp.KSubList(np.ones(2), (1, 1), [0, 2], score=0.7)

        queue.push(sublist1)
        queue.push(sublist2)

        coords_1 = puw.quantity(np.array([
            [[1., 1., 1.],
             [2., 2., 2.],
             ]]), "angstroms")
        coords_2 = puw.quantity(np.array([
            [[5., 5., 5.],
             [6., 6., 6.],
             [7., 7., 7.],
             ],
            [[9., 9., 9.],
             [10., 10., 10.],
             [11., 11., 11.],
             ],
        ]), "angstroms")

        feature_lists = [
            cp.FeatureList("AD", distances=np.ones(2), coords=coords_1),
            cp.FeatureList("AAR", distances=np.ones(4), coords=coords_2),
        ]

        pharmacophores = cp.CommonPharmacophoreFinder._get_pharmacophores(queue, feature_lists)
        radius = puw.quantity(1.0, "angstroms")
        expected_1 = Pharmacophore([
            PharmacophoricPoint("hb acceptor", coords_2[1][0], radius),
            PharmacophoricPoint("aromatic ring", coords_2[1][2], radius),
        ], ref_mol=1, ref_struct=1, score=0.7)
        expected_2 = Pharmacophore([
            PharmacophoricPoint("hb acceptor", coords_1[0][0], radius),
            PharmacophoricPoint("hb donor", coords_1[0][1], radius),
        ], ref_mol=0, ref_struct=0, score=0.5)
        assert pharmacophores == [expected_1, expected_2]


class TestFeatList:

    def test_distance_vector(self):
        coords = np.array([
            [[4., 3., 0.],
             [4., 0., 0.],
             [0., 0., 0.],
             ],
            [[3.125, 3.903123749, 0.],
             [4., 0., 0.],
             [0., 0., 0.],
             ],
        ])
        distance_vector = cp.FeatureList.distance_vector(coords)
        expected = np.array([
            [3., 5., 4.],
            [4., 5., 4.],
        ])
        assert np.allclose(distance_vector, expected)

    def test_has_variant(self):
        feat_list = cp.FeatureList(variant="AADPR", distances=None, coords=None)
        assert feat_list.has_variant("AAR")
        assert not feat_list.has_variant("AHR")

    def test_k_sublists(self):
        coords = puw.quantity(np.ones((2, 5, 3)), "angstroms")
        distances = np.array([
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
        ])
        feat_list = cp.FeatureList(variant="AADPR", distances=distances, coords=coords)
        k_variant = cp.KVariant((0, 2, 3), 1)
        k_sublists = feat_list.k_sublists(k_variant)
        assert k_sublists == [
            cp.KSubList(np.array([1, 2, 7]), (1, 0), [0, 2, 3]),
            cp.KSubList(np.array([11, 12, 17]), (1, 1), [0, 2, 3]),
        ]

    def test_sub_variant_distances(self):
        distances = np.array([
            [0, 1, 2, 3, 4,
             5, 6, 7, 8,
             9, 10, 11,
             12, 13,
             14
             ],
        ])
        feat_list = cp.FeatureList(variant="ABCDEF", distances=distances, coords=None)
        sub_variant = [1, 3, 4, 5]  # Subvariant BDEF
        distances = feat_list.subvariant_distances(sub_variant, conf=0)
        expected = np.array([
            6, 7, 8,
            12, 13,
            14,
        ])
        assert np.all(distances == expected)


class TestKSubList:

    def test_id_is_incremented(self):
        cp.KSubList.ID = 1

        distance = np.array([1, 2, 30])
        sublist_1 = cp.KSubList(distance, (0, 0), [1, 2, 3, 4])
        sublist_2 = cp.KSubList(distance, (0, 0), [1, 2, 3, 4])
        sublist_3 = cp.KSubList(distance, (0, 0), [1, 2, 3, 4])

        assert sublist_1.id == 1
        assert sublist_2.id == 2
        assert sublist_3.id == 3


class TestSurvivingBox:

    def test_equality(self):
        box_1 = cp.SurvivingBox()
        box_1.id = [1, 2, 3]

        box_2 = cp.SurvivingBox()
        box_2.id = [2, 3]

        assert not box_1 == box_2

        box_3 = cp.SurvivingBox()
        box_3.id = [1, 2, 3]
        assert box_1 == box_3


class TestScoringFunction:

    def test_point_score(self):
        ref_coords = puw.quantity(np.array([
            [4., 3., 0.],
            [4., 0., 0.],
            [0., 0., 0.]
        ]), "angstroms")
        coords = puw.quantity(np.array([
            [6.3333333333, 2.98142397, 0.],
            [6., 0., 0.],
            [0., 0., 0.],
        ]), "angstroms")

        scoring_fn = cp.ScoringFunction()
        score = scoring_fn.point_score(align_pharmacophores(ref_coords, coords))
        assert score == 0.14632288133597426


class TestFeatListQueue:

    def test_gives_priority_to_higher_scores(self):
        queue = cp.FeatListQueue(size=2)
        sublist1 = cp.KSubList(np.array([1, 1, 1]), (0, 0), [0, 2, 3], score=0.5)
        sublist2 = cp.KSubList(np.array([2, 2, 2]), (1, 0), [0, 2, 3], score=0.7)
        sublist3 = cp.KSubList(np.array([3, 3, 3]), (2, 0), [0, 2, 3], score=0.8)

        queue.push(sublist1)
        queue.push(sublist2)
        queue.push(sublist3)

        assert queue.pop() == sublist2
        assert queue.pop() == sublist3
        assert queue.empty()

    def test_cannot_add_items_with_lower_score_to_full_queue(self):
        queue = cp.FeatListQueue(size=2)
        sublist1 = cp.KSubList(np.array([1, 1, 1]), (0, 0), [0, 2, 3], score=0.5)
        sublist2 = cp.KSubList(np.array([2, 2, 2]), (1, 0), [0, 2, 3], score=0.7)
        sublist3 = cp.KSubList(np.array([3, 3, 3]), (2, 0), [0, 2, 3], score=0.3)

        queue.push(sublist1)
        queue.push(sublist2)
        queue.push(sublist3)  # Should not be pushed

        assert queue.pop() == sublist1
        assert queue.pop() == sublist2
        assert queue.empty()

    def test_peek(self):
        queue = cp.FeatListQueue(size=7)
        sublist2 = cp.KSubList(np.array([2, 2, 2]), (1, 0), [0, 2, 3], score=0.8)
        sublist1 = cp.KSubList(np.array([1, 1, 1]), (0, 0), [0, 2, 3], score=0.5)
        sublist3 = cp.KSubList(np.array([3, 3, 3]), (2, 0), [0, 2, 3], score=0.9)

        queue.push(sublist1)
        queue.push(sublist2)
        queue.push(sublist3)

        assert queue.peek_left().score == 0.5
        assert queue.peek_right().score == 0.9

    def test_cannot_add_repeated_elements(self):
        queue = cp.FeatListQueue(size=2)
        sublist1 = cp.KSubList(np.array([1, 1, 1]), (0, 0), [0, 2, 3], score=0.5)
        sublist2 = cp.KSubList(np.array([2, 2, 2]), (1, 0), [0, 2, 3], score=0.7)

        queue.push(sublist1)
        queue.push(sublist2)
        queue.push(sublist2)

        assert queue.pop() == sublist1
        assert queue.pop() == sublist2
        assert queue.empty()
