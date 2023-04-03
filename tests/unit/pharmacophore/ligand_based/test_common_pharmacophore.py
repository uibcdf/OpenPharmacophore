import pyunitwizard as puw
import numpy as np
from openpharmacophore import Pharmacophore, PharmacophoricPoint
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

        expected_bins = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        assert finder.min_dist == min_dist
        assert finder.max_dist == max_dist
        assert finder.bin_size == bin_size
        assert np.all(finder.bins == expected_bins)

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
        cp_finder = cp.CommonPharmacophoreFinder(n_points=3)
        common_pharmacophores = cp_finder(molecules)

        radius = puw.quantity(1.0, "angstroms")
        assert common_pharmacophores == [
            Pharmacophore([
                PharmacophoricPoint("hb acceptor", coords_1[0], radius),
                PharmacophoricPoint("hb donor", coords_1[0], radius),
                PharmacophoricPoint("hydrophobicity", coords_1[0], radius),
            ],
                ref_struct=0,
                ref_mol=0
            ),
        ]

    def test_get_feat_lists(self):
        coords_1 = puw.quantity(np.array([
            [4., 0., 0.],
            [0., 0., 0.],
        ]), "angstroms")

        conformer_1 = ChemFeatContainer([
            ChemFeat(coords_1[0], "hb acceptor"),
            ChemFeat(coords_1[1], "hb donor"),
        ])

        coords_2 = puw.quantity(np.array([
            [6., 0., 0.],
            [0., 0., 0.],
        ]), "angstroms")

        conformer_2 = ChemFeatContainer([
            ChemFeat(coords_2[0], "aromatic ring"),
            ChemFeat(coords_2[1], "hb donor"),
        ])

        cp_finder = cp.CommonPharmacophoreFinder(n_points=3)
        feat_lists = cp_finder._get_feat_lists([
            [conformer_1], [conformer_2]
        ])
        assert len(feat_lists) == 2
        assert feat_lists[0].variant == "AD"
        assert feat_lists[0].mol_id == (0, 0)
        assert np.all(feat_lists[0].distances == np.array([4.]))
        assert np.all(feat_lists[0].coords == coords_1)

        assert feat_lists[1].variant == "DR"
        assert feat_lists[1].mol_id == (1, 0)
        assert np.all(feat_lists[1].distances == np.array([6.]))
        assert np.all(feat_lists[1].coords == coords_2)

    def test_common_k_point_variants(self):
        coords = puw.quantity([1., 2., 3.], "angstroms")
        distances = np.array([1., 2, .3])
        all_lists = [
            cp.FeatureList("ADHP", (0, 0), distances, coords),
            cp.FeatureList("ADHR", (1, 0), distances, coords),
            cp.FeatureList("ADHN", (1, 1), distances, coords),
        ]


class TestFeatList:

    def test_distance_vector(self):
        coords = np.array([
            [4., 3., 0.],
            [4., 0., 0.],
            [0., 0., 0.],
        ])
        distance_vector = cp.FeatureList.distance_vector(coords)
        assert np.all(distance_vector == np.array([3., 5., 4.]))
