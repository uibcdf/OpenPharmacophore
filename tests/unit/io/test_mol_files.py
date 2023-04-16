from fake_files import fake_mol2_file
from openpharmacophore.io.mol_files import _mol2, _iter_mol2


class TestMol2:

    def test_load_mol2_file(self):
        molecules = _mol2(fake_mol2_file())

        assert len(molecules) == 2
        assert molecules[0].GetNumAtoms() == 14
        assert molecules[1].GetNumAtoms() == 25

    def test_mol2_iterator(self):
        molecules = []
        for mol in _iter_mol2(fake_mol2_file()):
            molecules.append(mol)

        assert len(molecules) == 2
        assert molecules[0].GetNumAtoms() == 14
        assert molecules[1].GetNumAtoms() == 25
