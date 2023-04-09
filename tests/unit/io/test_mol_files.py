from openpharmacophore.io.mol_files import mol2, iter_mol2


class TestMol2:

    def test_load_mol2_file(self, ligands_mol2):
        molecules = mol2(file_name=ligands_mol2)

        assert len(molecules) == 3
        assert molecules[0].GetNumAtoms() == 14
        assert molecules[1].GetNumAtoms() == 25
        assert molecules[2].GetNumAtoms() == 29

    def test_mol2_iterator(self, ligands_mol2):
        molecules = []
        for mol in iter_mol2(ligands_mol2):
            molecules.append(mol)

        assert len(molecules) == 3
        assert molecules[0].GetNumAtoms() == 14
        assert molecules[1].GetNumAtoms() == 25
        assert molecules[2].GetNumAtoms() == 29
