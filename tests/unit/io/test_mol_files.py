from io import StringIO
from rdkit import Chem
from fake_files import fake_mol2_file, fake_sdf_file
from openpharmacophore import mol_files


class TestMol2:

    def test_load_mol2_file(self):
        molecules = mol_files._mol2(fake_mol2_file())

        assert len(molecules) == 2
        assert molecules[0].GetNumAtoms() == 14
        assert molecules[1].GetNumAtoms() == 25

    def test_mol2_iterator(self):
        molecules = []
        for mol in mol_files._iter_mol2(fake_mol2_file()):
            molecules.append(mol)

        assert len(molecules) == 2
        assert molecules[0].GetNumAtoms() == 14
        assert molecules[1].GetNumAtoms() == 25


class TestSmi:

    def test_file_with_header(self):
        contents = "This is the header\n" \
                   "CN=C=O\n" \
                   "COc1cc(C=O)ccc1O\n"
        file_ = StringIO(contents)
        molecules = mol_files._smi(file_, header=True)
        assert len(molecules) == 2
        assert molecules[0].GetNumAtoms() == 4
        assert molecules[1].GetNumAtoms() == 11

    def test_reads_mol_names(self):
        contents = "CN=C=O MIC\n" \
                   "COc1cc(C=O)ccc1O Vanillin\n"
        file_ = StringIO(contents)
        molecules = mol_files._smi(file_, header=False)
        assert len(molecules) == 2
        assert molecules[0].GetNumAtoms() == 4
        assert molecules[0].GetProp("_Name") == "MIC"
        assert molecules[1].GetNumAtoms() == 11
        assert molecules[1].GetProp("_Name") == "Vanillin"

    def test_iter_smi(self):
        contents = "This is the header\n" \
                   "CN=C=O MIC\n" \
                   "COc1cc(C=O)ccc1O Vanillin\n"
        file_ = StringIO(contents)
        molecules = [m for m in mol_files._iter_smi(file_, header=True)]
        assert len(molecules) == 2
        assert molecules[0].GetNumAtoms() == 4
        assert molecules[0].GetProp("_Name") == "MIC"
        assert molecules[1].GetNumAtoms() == 11
        assert molecules[1].GetProp("_Name") == "Vanillin"


class TestSdf:

    def test_read_sdf(self):
        molecules = mol_files._sdf(fake_sdf_file(), Chem.ForwardSDMolSupplier)
        assert len(molecules) == 2
        assert molecules[0].GetNumAtoms() == 75
        assert molecules[1].GetNumAtoms() == 62

    def test_iter_sdf(self):
        molecules = [
            m for m in mol_files._iter_sdf(fake_sdf_file(), Chem.ForwardSDMolSupplier)
        ]
        assert len(molecules) == 2
        assert molecules[0].GetNumAtoms() == 75
        assert molecules[1].GetNumAtoms() == 62
