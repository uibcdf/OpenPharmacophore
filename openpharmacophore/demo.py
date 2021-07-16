"""Module with demonstration objects and files.

This module contains useful objects and files to help document and demonstrate the use of
OpenPharmacophore.

Attributes
----------
pharmacophore: :obj:`openpharmacophore.Pharmacophore`
    Demonstration native pharmacophore object.

pharmacophore_pharmer_file : str
    Filename of a json file of a pharmacophore produced with `Pharmer <https://sourceforge.net/projects/pharmer/>`_.


Examples
-------
>>> import openpharmacophore as oph
>>> pharmacophore = oph.demo.pharmacophore


Todo
----
These demonstration objects and files should be replaced. The current ones are temporary:
    - The pharmer file is coming from a third repository. We should produce our own file.
    - The native pharmacophore is built on the fly in this module. But should be just read from a
      native file.

"""


import pkg_resources
from openpharmacophore import pharmacophoric_elements as elements
from openpharmacophore import Pharmacophore

pharmacophore_pharmer_file = pkg_resources.resource_filename('openpharmacophore', 'data/pharmer.json')

pharmacophore = Pharmacophore()
#pharmacophore.add_element(elements.PositiveChargeSphere('[0,0,0] angstroms', '1.0 angstroms'))
#pharmacophore.add_element(elements.NegativeChargeSphere('[-1,2,0] angstroms', '1.0 angstroms'))
#pharmacophore.add_element(elements.HBAcceptorSphereAndVector('[-1,-1,0] angstroms', '1.0 angstroms',[-1,-1,-2]))
#pharmacophore.add_element(elements.HydrophobicGaussianKernel('[1,1,3] angstroms', '1.0 angstroms'))
#pharmacophore.add_element(elements.AromaticRingSphere('[-2,-3,0] angstroms', '1.5 angstroms'))

del(pkg_resources, elements, Pharmacophore)

