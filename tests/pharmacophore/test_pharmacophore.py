from openpharmacophore.pharmacophore.pharmacophore import Pharmacophore
from openpharmacophore import LigandBasedPharmacophore, PharmacophoricPoint
import pytest
import pyunitwizard as puw


def test_cannot_instantiate_pharmacophore():

    with pytest.raises(TypeError):
        ph = Pharmacophore()
