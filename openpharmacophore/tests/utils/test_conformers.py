import openpharmacophore.data as data
from openpharmacophore.utils.conformers import generate_conformers, conformer_energy
from rdkit import Chem


def test_generate_conformers():

    supplier = Chem.SmilesMolSupplier(data.ligands["mols"])
    n_conformers = [1, 2, 3, 2, 1]
    for ii, mol in enumerate(supplier):
        mol = generate_conformers(mol, n_conformers[ii])
        assert mol.GetNumConformers() == n_conformers[ii]


def test_conformers_energy():

    supplier = Chem.SmilesMolSupplier(data.ligands["mols"])
    molecule = supplier[0]
    molecule = generate_conformers(molecule, 2)

    energy_conformer_1 = conformer_energy(molecule, conformer_id=0,
                                          forcefield="UFF")
    energy_conformer_2 = conformer_energy(molecule, conformer_id=1,
                                          forcefield="MMFF")
    assert isinstance(energy_conformer_1, float)
    assert isinstance(energy_conformer_2, float)
    assert energy_conformer_1 != energy_conformer_2
