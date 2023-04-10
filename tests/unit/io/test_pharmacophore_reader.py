from openpharmacophore import pharmacophore_reader
import numpy as np
import pyunitwizard as puw


class TestReadMol2:

    def test_read_mol2(self, mol2_pharmacophore_path):
        pharmacophores = pharmacophore_reader.read_mol2(mol2_pharmacophore_path)
        assert len(pharmacophores) == 2

        pharma_1 = pharmacophores[0]
        expected_feats = ["hydrophobicity", "positive charge", "aromatic ring"]
        expected_centers = [
            np.array([5.3180, 0.5997, 0.6328]),
            np.array([7.9404, -0.6298, 0.7958]),
            np.array([-1.5243, -0.3433, -1.1795]),
        ]
        expected_radius = [1.0, 1.5, 1.0]

        assert [p.feature_name for p in pharma_1] == expected_feats
        for ii in range(len(pharma_1)):
            assert np.allclose(puw.get_value(pharma_1[ii].center), expected_centers[ii])
            assert np.allclose(puw.get_value(pharma_1[ii].radius), expected_radius[ii])

        assert pharma_1.score == 0.8452
        assert pharma_1.ref_mol == 1
        assert pharma_1.ref_struct == 2
        assert pharma_1.props["min_actives"] == "4"

        pharma_2 = pharmacophores[1]
        expected_feats = ["hb acceptor", "hydrophobicity", "positive charge"]
        expected_centers = [
            np.array([-0.4548, 0.0480, 0.2469]),
            np.array([-4.7914, -1.2623, 1.4397]),
            np.array([-7.1862, -2.5267, 0.5723]),
        ]
        expected_radius = [1.0, 2.0, 1.45]

        assert [p.feature_name for p in pharma_2] == expected_feats
        for ii in range(len(pharma_2)):
            assert np.allclose(puw.get_value(pharma_2[ii].center), expected_centers[ii])
            assert np.allclose(puw.get_value(pharma_2[ii].radius), expected_radius[ii])

        assert pharma_2.score is None
        assert pharma_2.ref_mol is None
        assert pharma_2.ref_struct is None
        assert len(pharma_2.props) == 0


class TestReadPH4:

    def test_read_ph4(
            self,
            moe_pharmacophore_path
    ):
        points = pharmacophore_reader.read_ph4(moe_pharmacophore_path)
        points = sorted(points, key=lambda p: p.short_name)
        assert len(points) == 10

        # HB Acceptor features should be first
        acceptor = points[0]
        assert acceptor.feature_name == "hb acceptor"
        assert np.all(
            np.around(puw.get_value(acceptor.center, "angstroms"), 2) == np.around(np.array([0.312, 3.0175, -2.44825]),
                                                                                   2)
        )
        assert np.around(puw.get_value(acceptor.radius, "angstroms"), 2) == 0.57
        acceptor = points[1]
        assert acceptor.feature_name == "hb acceptor"
        assert np.all(
            np.around(puw.get_value(acceptor.center, "angstroms"), 2) ==
            np.around(np.array([-1.95875, 2.536, -3.03625]), 2)
        )
        assert np.around(puw.get_value(acceptor.radius, "angstroms"), 2) == 0.62
        acceptor_2 = points[2]
        assert acceptor_2.feature_name == "hb acceptor"
        assert np.all(
            np.around(puw.get_value(acceptor_2.center, "angstroms"), 2) ==
            np.around(np.array([-0.755095833333333, 6.3286375, -3.96758333333333]), 2)
        )
        assert np.around(puw.get_value(acceptor_2.radius, "angstroms"), 2) == 1.25

        # Donor points
        donor = points[3]
        assert donor.feature_name == "hb donor"
        assert np.all(
            np.around(puw.get_value(donor.center, "angstroms"), 2) == np.around(np.array([1.71, 1.43075, -1.4255]), 2)
        )
        assert np.around(puw.get_value(donor.radius, "angstroms"), 2) == 0.51

        # Hydrophobic points
        hyd = points[4]
        assert hyd.feature_name == "hydrophobicity"
        assert np.all(
            np.around(puw.get_value(hyd.center, "angstroms"), 2) == np.around(np.array([2.7895, 2.4035, -1.40875]), 2)
        )
        assert np.around(puw.get_value(hyd.radius, "angstroms"), 2) == 0.55
        hyd_2 = points[5]
        assert hyd_2.feature_name == "hydrophobicity"
        assert np.all(
            np.around(puw.get_value(hyd_2.center, "angstroms"), 2) == np.around(
                np.array([-1.54725, -2.979375, -0.961875]),
                2)
        )
        assert np.around(puw.get_value(hyd_2.radius, "angstroms"), 2) == 0.74

        # Aromatic Points
        aromatic_1 = points[6]
        assert aromatic_1.feature_name == "aromatic ring"
        assert np.all(
            np.around(puw.get_value(aromatic_1.center, "angstroms"), 2) ==
            np.around(np.array([-0.748458333333333, 2.13108333333333, -2.490375]), 2)
        )
        assert np.around(puw.get_value(aromatic_1.radius, "angstroms"), 2) == 0.58
        aromatic_2 = points[7]
        assert aromatic_2.feature_name == "aromatic ring"
        assert np.all(
            np.around(puw.get_value(aromatic_2.center, "angstroms"), 2) ==
            np.around(np.array([-1.719625, -0.0273333333333334, -2.055625]), 2)
        )
        assert np.around(puw.get_value(aromatic_2.radius, "angstroms"), 2) == 0.6
        aromatic_3 = points[8]
        assert aromatic_3.feature_name == "aromatic ring"
        assert np.all(
            np.around(puw.get_value(aromatic_3.center, "angstroms"), 2) ==
            np.around(np.array([5.20029166666667, 1.25479166666667, -0.199041666666667]), 2)
        )
        assert np.around(puw.get_value(aromatic_3.radius, "angstroms"), 2) == 0.61
        aromatic_4 = points[9]
        assert aromatic_4.feature_name == "aromatic ring"
        assert np.all(
            np.around(puw.get_value(aromatic_4.center, "angstroms"), 2) ==
            np.around(np.array([-0.755095833333333, 6.3286375, -3.96758333333333]), 2)
        )
        assert np.around(puw.get_value(aromatic_4.radius, "angstroms"), 2) == 1.25


class TestReadJSON:

    def test_read_json(self, json_pharmacophore_path):
        pharmacophore = pharmacophore_reader.read_json(json_pharmacophore_path)

        assert len(pharmacophore) == 5

        assert pharmacophore[0].feature_name == "hb acceptor"
        assert pharmacophore[0].has_direction
        assert puw.get_value(pharmacophore[0].radius, "angstroms") == 1.0
        assert np.allclose(puw.get_value(pharmacophore[0].center, "angstroms"),
                           np.array([21.352, -14.531, 19.625]))
        assert np.allclose(pharmacophore[0].direction,
                           np.array([-0.6405836470264256, 0.7029084735090229, -0.3091476492414897]))

        assert pharmacophore[1].feature_name == "hb acceptor"
        assert pharmacophore[1].has_direction
        assert puw.get_value(pharmacophore[1].radius, "angstroms") == 1.0
        assert np.allclose(puw.get_value(pharmacophore[1].center, "angstroms"),
                           np.array([19.355, -18.32, 23.987]))
        assert np.allclose(pharmacophore[1].direction,
                           np.array([0.6859059711903811, 0.09092493673854565, 0.721987295292979]))

        assert pharmacophore[2].feature_name == "hb donor"
        assert pharmacophore[2].has_direction
        assert puw.get_value(pharmacophore[2].radius, "angstroms") == 1.0
        assert np.allclose(puw.get_value(pharmacophore[2].center, "angstroms"),
                           np.array([20.977, -16.951, 18.746]))
        assert np.allclose(pharmacophore[2].direction,
                           np.array([0.71662539652105, -0.5202950182802607, -0.46447942367105033]))

        assert pharmacophore[3].feature_name == "negative charge"
        assert not pharmacophore[3].has_direction
        assert puw.get_value(pharmacophore[3].radius, "angstroms") == 1.5
        assert np.allclose(puw.get_value(pharmacophore[3].center, "angstroms"),
                           np.array([21.66899, -15.077667, 20.608334]))

        assert pharmacophore[4].feature_name == "negative charge"
        assert not pharmacophore[4].has_direction
        assert puw.get_value(pharmacophore[4].radius, "angstroms") == 2.0
        assert np.allclose(puw.get_value(pharmacophore[4].center, "angstroms"),
                           np.array([19.985, -19.404402, 22.8422]))


class TestReadLigandScout:

    def test_from_ligandscout(
            self,
            ligand_scout_pharmacophore_path
    ):
        points = pharmacophore_reader.read_ligandscout(ligand_scout_pharmacophore_path)
        assert len(points) == 4

        neg_ion = points[0]
        assert neg_ion.feature_name == "negative charge"
        assert np.all(
            np.around(puw.get_value(neg_ion.center, "angstroms"), 1) == np.array([-8.0, 10.0, -9.5])
        )
        assert puw.get_value(neg_ion.radius, "angstroms") == 1.5

        donor = points[1]
        assert donor.feature_name == "hb donor"
        assert np.all(
            np.around(puw.get_value(donor.center, "angstroms"), 1) == np.array([-8.0, 2.0, -10.0])
        )
        assert puw.get_value(donor.radius, "angstroms") == 1.5
        dir_expected = np.array([-5.690445899963379, 0.5822541117668152, -10.5515718460083])
        dir_expected /= np.linalg.norm(dir_expected)
        assert np.all(
            np.around(donor.direction, 2) == np.around(dir_expected, 2)
        )

        ring = points[2]
        assert ring.feature_name == "aromatic ring"
        assert np.all(
            np.around(puw.get_value(ring.center, "angstroms"), 1) == np.array([0.0, 6.5, -3.0])
        )
        assert puw.get_value(ring.radius, "angstroms") == 1.5
        dir_expected = np.array([-3.8126893043518066, 1.7578959465026855, 0.6093783378601074])
        dir_expected /= np.linalg.norm(dir_expected)
        assert np.all(
            np.around(ring.direction, 2) == np.around(dir_expected, 2)
        )

        excluded_vol = points[3]
        assert excluded_vol.feature_name == "excluded volume"
        assert np.all(
            np.around(puw.get_value(excluded_vol.center, "angstroms"), 1) == np.array([5.5, 4.5, -2.0])
        )
        assert round(puw.get_value(excluded_vol.radius, "angstroms"), 1) == 1.0
