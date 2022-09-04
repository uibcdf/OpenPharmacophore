from openpharmacophore.io.pharmacophore_mol2 import (
    load_mol2_pharmacophoric_points, mol2_file_info, pad_coordinate_with_zeros)
import openpharmacophore.data as data
from example_pharmacophores import two_element_pharmacophore, five_element_pharmacophore


def test_load_mol2_pharmacophoric_points():
    streptadivin_file = data.pharmacophores["streptadivin"]
    pharmacophores_ = load_mol2_pharmacophoric_points(streptadivin_file)
    assert len(pharmacophores_) == 6
    assert len(pharmacophores_[0]) == 9
    assert len(pharmacophores_[1]) == 10
    assert len(pharmacophores_[2]) == 16
    assert len(pharmacophores_[3]) == 11
    assert len(pharmacophores_[4]) == 13
    assert len(pharmacophores_[5]) == 13

    elastase_file = data.pharmacophores["elastase"]
    pharmacophores_ = load_mol2_pharmacophoric_points(elastase_file)

    assert len(pharmacophores_) == 8
    assert len(pharmacophores_[0]) == 4
    assert len(pharmacophores_[1]) == 21
    assert len(pharmacophores_[2]) == 16
    assert len(pharmacophores_[3]) == 18
    assert len(pharmacophores_[4]) == 10
    assert len(pharmacophores_[5]) == 12
    assert len(pharmacophores_[6]) == 11
    assert len(pharmacophores_[7]) == 14


def test_pad_coordinate_with_zeros_positive_number():
    assert pad_coordinate_with_zeros(25) == "25.000"
    assert pad_coordinate_with_zeros(2.5) == "2.5000"
    assert pad_coordinate_with_zeros(0.25) == "0.2500"


def test_pad_coordinate_with_zeros_negative_number():
    assert pad_coordinate_with_zeros(-25) == "-25.000"
    assert pad_coordinate_with_zeros(-2.5) == "-2.5000"
    assert pad_coordinate_with_zeros(-0.25) == "-0.2500"


def test_mol2_file_info():
    # Test for two element pharmacophore
    mol2_list = mol2_file_info(two_element_pharmacophore())
    expected_output = ['@<TRIPOS>MOLECULE\n',
                       '@<TRIPOS>ATOM\n',
                       '      1 AR           1.0000    0.0000    0.0000   AR     0   AR      0.0000\n',
                       '      2 ACC          1.0000    2.0000    2.0000   HB     1   HB      0.0000\n',
                       '@<TRIPOS>BOND\n']
    assert mol2_list == expected_output

    # Test for five element pharmacophore
    mol2_list = mol2_file_info(five_element_pharmacophore())
    expected_output = ['@<TRIPOS>MOLECULE\n',
                       '@<TRIPOS>ATOM\n',
                       '      1 ACC          1.0000    2.0000    2.0000   HB     0   HB      0.0000\n',
                       '      2 DON          1.0000    2.0000    2.0000   HB     1   HB      0.0000\n',
                       '      3 HYD         -1.0000    2.0000    2.0000   HYD    2   HYD     0.0000\n',
                       '      4 AR           1.0000    0.0000    0.0000   AR     3   AR      0.0000\n',
                       '      5 AR           0.0000    1.0000    2.0000   AR     4   AR      0.0000\n',
                       '@<TRIPOS>BOND\n']
    assert mol2_list == expected_output
