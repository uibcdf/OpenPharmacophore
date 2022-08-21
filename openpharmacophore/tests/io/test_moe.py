import openpharmacophore.io as io
import openpharmacophore.data as data
import numpy as np
import pyunitwizard as puw
import datetime
from example_pharmacophores import three_element_pharmacophore


def test_from_moe():
    file_name = data.pharmacophores["gmp"]
    points = io.from_moe(file_name)
    assert len(points) == 10

    # HB Acceptor features should be first
    acceptor = points[0]
    assert acceptor.feature_name == "hb acceptor"
    assert np.all(
        np.around(puw.get_value(acceptor.center, "angstroms"), 2) == np.around(np.array([0.312, 3.0175, -2.44825]), 2)
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
        np.around(puw.get_value(hyd_2.center, "angstroms"), 2) == np.around(np.array([-1.54725, -2.979375, -0.961875]),
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


def test_to_moe():
    pharmacophore_str = io._moe_ph4_string(three_element_pharmacophore())

    now = datetime.datetime.now()
    month = str(now.month)
    year = str(now.year)
    expected_str = (f'#moe:ph4que {year}.{month}\n'
                    f'#pharmacophore 5 tag t value *\n'
                    f'scheme t Unified matchsize i 0 title t s $\n'
                    f'#feature 3 expr tt color ix x r y r z r r r ebits ix gbits ix\n'
                    f'Acc df2f2 0.9999999999999999 1.9999999999999998 1.9999999999999998 0.9999999999999999 0 300 Aro df2f2 0.9999999999999999 0.0 0.0 0.9999999999999999 0 300 \n'
                    f'#volumesphere 90 x r y r z r r r\n'
                    f'1.9999999999999998 0.9999999999999999 1.9999999999999998 0.9999999999999999 \n'
                    f'#endpharmacophore')

    assert pharmacophore_str == expected_str
