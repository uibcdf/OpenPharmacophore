from openpharmacophore import pharmacophore_reader
import numpy as np
import pyunitwizard as puw


class TestReadMol2:

    def test_read_mol2(
            self,
            mol2_pharmacophore_path_elastase,
            mol2_pharmacophore_path_streptadivin
    ):
        pharmacophores_ = pharmacophore_reader.read_mol2(mol2_pharmacophore_path_streptadivin)
        # TODO: we can test with a smaller file
        assert len(pharmacophores_) == 6
        assert len(pharmacophores_[0]) == 9
        assert len(pharmacophores_[1]) == 10
        assert len(pharmacophores_[2]) == 16
        assert len(pharmacophores_[3]) == 11
        assert len(pharmacophores_[4]) == 13
        assert len(pharmacophores_[5]) == 13

        pharmacophores_ = pharmacophore_reader.read_mol2(mol2_pharmacophore_path_elastase)

        assert len(pharmacophores_) == 8
        assert len(pharmacophores_[0]) == 4
        assert len(pharmacophores_[1]) == 21
        assert len(pharmacophores_[2]) == 16
        assert len(pharmacophores_[3]) == 18
        assert len(pharmacophores_[4]) == 10
        assert len(pharmacophores_[5]) == 12
        assert len(pharmacophores_[6]) == 11
        assert len(pharmacophores_[7]) == 14


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
        points, molecular_system, ligand = pharmacophore_reader.read_json(
            json_pharmacophore_path, load_mol_sys=False)

        assert len(points) == 5
        assert molecular_system is None
        assert ligand is None

        assert points[0].feature_name == "hb acceptor"
        assert points[0].has_direction
        assert puw.get_value(points[0].radius, "angstroms") == 1.0
        assert np.allclose(puw.get_value(points[0].center, "angstroms"),
                           np.array([21.352, -14.531, 19.625]))
        assert np.allclose(points[0].direction,
                           np.array([-0.6405836470264256, 0.7029084735090229, -0.3091476492414897]))

        assert points[1].feature_name == "hb acceptor"
        assert points[1].has_direction
        assert puw.get_value(points[1].radius, "angstroms") == 1.0
        assert np.allclose(puw.get_value(points[1].center, "angstroms"),
                           np.array([19.355, -18.32, 23.987]))
        assert np.allclose(points[1].direction,
                           np.array([0.6859059711903811, 0.09092493673854565, 0.721987295292979]))

        assert points[2].feature_name == "hb donor"
        assert points[2].has_direction
        assert puw.get_value(points[2].radius, "angstroms") == 1.0
        assert np.allclose(puw.get_value(points[2].center, "angstroms"),
                           np.array([20.977, -16.951, 18.746]))
        assert np.allclose(points[2].direction,
                           np.array([0.71662539652105, -0.5202950182802607, -0.46447942367105033]))

        assert points[3].feature_name == "negative charge"
        assert not points[3].has_direction
        assert puw.get_value(points[3].radius, "angstroms") == 1.5
        assert np.allclose(puw.get_value(points[3].center, "angstroms"),
                           np.array([21.66899, -15.077667, 20.608334]))

        assert points[4].feature_name == "negative charge"
        assert not points[4].has_direction
        assert puw.get_value(points[4].radius, "angstroms") == 2.0
        assert np.allclose(puw.get_value(points[4].center, "angstroms"),
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
