from io import StringIO
from rdkit import Chem
from openpharmacophore import mol_db


class TestInMemoryDB:

    def test_iterate_smiles_list(self):
        smiles = [
            "CN=C=O",
            "COc1cc(C=O)ccc1O",
        ]
        db = mol_db.InMemoryMolDB()
        db.from_smiles(smiles)

        molecules = []
        for mol in db:
            molecules.append(mol)

        assert len(molecules) == 2
        assert molecules[0].n_atoms == 4
        assert molecules[1].n_atoms == 11

    def test_iterate_rdkit_mol(self):
        rdkit_mols = [
            Chem.MolFromSmiles("CN=C=O"),
            Chem.MolFromSmiles("COc1cc(C=O)ccc1O"),
        ]

        db = mol_db.InMemoryMolDB()
        db.from_rdkit(rdkit_mols)

        molecules = []
        for mol in db:
            molecules.append(mol)

        assert len(molecules) == 2
        assert molecules[0].n_atoms == 4
        assert molecules[1].n_atoms == 11


class TestMolDB:

    def test_iterate_file_list(self):
        db = mol_db.MolDB()
        db._from_file(StringIO("CN=C=O\nCOc1cc(C=O)ccc1O"), "smi")
        db._from_file(StringIO("CN1CCC[C@H]1c2cccnc2"), "smi")

        molecules = []
        for mol in db:
            molecules.append(mol)

        assert len(molecules) == 3
        assert molecules[0].n_atoms == 4
        assert molecules[1].n_atoms == 11
        assert molecules[2].n_atoms == 12
