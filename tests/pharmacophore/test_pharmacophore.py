from openpharmacophore.pharmacophore.pharmacophore import Pharmacophore
import pytest


def test_cannot_instantiate_pharmacophore():

    with pytest.raises(TypeError):
        ph = Pharmacophore()
