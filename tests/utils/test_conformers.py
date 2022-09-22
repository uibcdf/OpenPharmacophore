from openpharmacophore.utils.conformers import generate_conformers
from rdkit.Chem import MolFromSmiles


def test_generate_conformers():

    molecule = MolFromSmiles("c1ccccc1")
    molecule = generate_conformers(molecule, 2)
    assert molecule.GetNumConformers() == 2
