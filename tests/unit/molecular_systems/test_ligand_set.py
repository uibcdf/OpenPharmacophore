import pytest
from openpharmacophore.molecular_systems.ligand import AbstractLigand, LigandSet


class FakeLigand(AbstractLigand):

    def __init__(self):
        self._hydrogens = False
        self._conformers = 0

    @property
    def has_hydrogens(self) -> bool:
        return self._hydrogens

    @property
    def n_conformers(self):
        return self._conformers

    def add_hydrogens(self):
        self._hydrogens = True

    def generate_conformers(self, n_confs):
        self._conformers = n_confs


def fake_ligand_set():
    ligs = LigandSet()
    for ii in range(3):
        ligs.add(FakeLigand())
    return ligs


def test_add_hydrogens_to_all_ligands():
    lig_set = fake_ligand_set()
    lig_set.add_hydrogens("all")

    assert all(lig.has_hydrogens for lig in lig_set)


def test_add_hydrogens_to_subset():
    lig_set = fake_ligand_set()
    lig_set.add_hydrogens([0, 2])

    assert lig_set[0].has_hydrogens
    assert not lig_set[1].has_hydrogens
    assert lig_set[2].has_hydrogens


def test_add_conformers_to_all_ligands():
    lig_set = fake_ligand_set()
    lig_set.generate_conformers(indices="all", n_confs=5)

    assert all([lig.n_conformers == 5 for lig in lig_set])


def test_add_conformers_to_selection():
    lig_set = fake_ligand_set()
    lig_set.generate_conformers(indices=[0, 1], n_confs=5)

    assert lig_set[0].n_conformers == 5
    assert lig_set[1].n_conformers == 5
    assert lig_set[2].n_conformers == 0
