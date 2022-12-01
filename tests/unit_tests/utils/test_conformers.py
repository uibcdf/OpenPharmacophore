from openpharmacophore.utils.conformers import generate_conformers
from rdkit.Chem import AddHs, MolFromSmiles


def test_generate_conformers_mol_with_no_hyd():
    molecule = MolFromSmiles("CN[C@H](Cc1ccccc1)C(=O)N2CCC[C@H]2C(=O)NCC3CCC(CC3)N")
    molecule = generate_conformers(molecule, 1)
    assert molecule.GetNumConformers() == 1
    assert any([a.GetSymbol() == "H" for a in molecule.GetAtoms()])


def test_generate_conformers_mol_with_hyd():
    molecule = MolFromSmiles("CN[C@H](Cc1ccccc1)C(=O)N2CCC[C@H]2C(=O)NCC3CCC(CC3)N")
    molecule = AddHs(molecule)
    assert molecule.GetNumConformers() == 0
    molecule = generate_conformers(molecule, 1)
    assert molecule.GetNumConformers() == 1
    assert any([a.GetSymbol() == "H" for a in molecule.GetAtoms()])
