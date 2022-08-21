from openpharmacophore import Pharmacophore
import openpharmacophore.io as io
import openpharmacophore.data as data
from example_pharmacophores import three_element_pharmacophore, five_element_pharmacophore


def test_read_pharmagist():
    streptadivin_file = data.pharmacophores["streptadivin"]
    pharmacophores_ = io.load_pharmacophores(streptadivin_file)
    pharmacophores_ = [Pharmacophore(ph) for ph in pharmacophores_]
    assert len(pharmacophores_) == 6
    for pharma in pharmacophores_:
        assert isinstance(pharma, Pharmacophore)

    assert len(pharmacophores_[0]) == 9
    assert len(pharmacophores_[1]) == 10
    assert len(pharmacophores_[2]) == 16
    assert len(pharmacophores_[3]) == 11
    assert len(pharmacophores_[4]) == 13
    assert len(pharmacophores_[5]) == 13

    elastase_file = data.pharmacophores["elastase"]
    pharmacophores_ = io.load_pharmacophores(elastase_file)
    pharmacophores_ = [Pharmacophore(ph) for ph in pharmacophores_]
    for pharma in pharmacophores_:
        assert isinstance(pharma, Pharmacophore)
    assert len(pharmacophores_) == 8
    assert len(pharmacophores_[0]) == 4
    assert len(pharmacophores_[1]) == 21
    assert len(pharmacophores_[2]) == 16
    assert len(pharmacophores_[3]) == 18
    assert len(pharmacophores_[4]) == 10
    assert len(pharmacophores_[5]) == 12
    assert len(pharmacophores_[6]) == 11
    assert len(pharmacophores_[7]) == 14


def test_to_pharmagist():
    # Test for two element pharmacophore
    mol2_list = io._pharmagist_file_info(three_element_pharmacophore())
    expected_output = ['@<TRIPOS>MOLECULE\n',
                       '@<TRIPOS>ATOM\n',
                       '      1 ACC          1.0000    2.0000    2.0000   HB     0   HB      0.0000\n',
                       '      2 AR           1.0000    0.0000    0.0000   AR     1   AR      0.0000\n',
                       '@<TRIPOS>BOND\n']
    assert mol2_list == expected_output

    # Test for five element pharmacophore
    mol2_list = io._pharmagist_file_info(five_element_pharmacophore())
    expected_output = ['@<TRIPOS>MOLECULE\n',
                       '@<TRIPOS>ATOM\n',
                       '      1 ACC          1.0000    2.0000    2.0000   HB     0   HB      0.0000\n',
                       '      2 DON          1.0000    2.0000    2.0000   HB     1   HB      0.0000\n',
                       '      3 HYD         -1.0000    2.0000    2.0000   HYD    2   HYD     0.0000\n',
                       '      4 AR           1.0000    0.0000    0.0000   AR     3   AR      0.0000\n',
                       '      5 AR           0.0000    1.0000    2.0000   AR     4   AR      0.0000\n',
                       '@<TRIPOS>BOND\n']
    assert mol2_list == expected_output
