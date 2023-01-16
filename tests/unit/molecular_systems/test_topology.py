from openpharmacophore.molecular_systems import Topology
from copy import deepcopy
import pytest
import numpy as np


@pytest.fixture
def topology_with_hydrogen():
    top = Topology()
    top.set_num_chains(2)
    top.add_atoms_to_chain(
        {
            "ALA": [("CA", "C"), ("O", "O"), ("NA", "N"), ("H", "H")],
            "MET": [("CA", "C"), ("O", "O"), ("CB", "C"), ("H", "H")],
        },
        0
    )
    top.add_atoms_to_chain(
        {
            "PRO": [("CA", "C"), ("O", "O"), ("NA", "N"), ("H", "H")],
            "LEU": [("CA", "C"), ("O", "O"), ("CB", "C"), ("H", "H")],
        },
        1
    )
    assert top.n_residues == 4
    assert top.n_atoms == 16
    return top


def test_create_topology_from_scratch():
    top = Topology()
    top.set_num_chains(2)
    assert top.n_chains == 2

    res_ind_1 = top.add_residue("ALA", chain=0)
    res_ind_2 = top.add_residue("MET", chain=1)
    assert res_ind_1 == 0
    assert res_ind_2 == 1
    assert top.n_residues == 2

    atom_ind_1 = top.add_atom("CA", "C", residue=res_ind_1)
    atom_ind_2 = top.add_atom("H", "H", residue=res_ind_2)
    assert atom_ind_1 == 0
    assert atom_ind_2 == 1
    assert top.n_atoms == 2


def test_has_hydrogens(topology_2_chains,
                       topology_with_hydrogen):
    assert not topology_2_chains.has_hydrogens()
    assert topology_with_hydrogen.has_hydrogens()


def test_has_ligands(
        topology_2_chains, topology_with_hydrogen,
        topology_with_ligand,
):
    assert not topology_2_chains.has_ligands()
    assert not topology_with_hydrogen.has_ligands()
    assert topology_with_ligand.has_ligands()


def test_returns_empty_list_when_topology_has_no_ligands(
        topology_2_chains
):
    lig_ids = topology_2_chains.ligand_ids()
    assert lig_ids == []


def test_find_all_ligands(topology_with_ligand):
    lig_ids = topology_with_ligand.ligand_ids()
    assert set(lig_ids) == {"EST:B", "EST:E"}


def test_get_residue_indices(topology_2_chains):
    expected_ind = [3, 4, 5]
    ind = topology_2_chains.get_residue_indices("MET", "A")
    assert ind == expected_ind


def test_get_residue_indices_extracts_correct_chain(topology_with_ligand):
    expected_ind = [16, 17, 18, 19]
    ind = topology_with_ligand.get_residue_indices("EST", "E")
    assert ind == expected_ind


def residue_in_topology(topology: Topology, res_name: str):
    for atom in topology.top.atoms:
        if atom.residue.name == res_name:
            return True
    return False


def test_remove_atoms(topology_with_ligand):
    new_topology = topology_with_ligand.remove_atoms(
        [8, 9, 10], inplace=False
    )
    assert new_topology.n_atoms == 17
    assert new_topology.n_chains == 4
    assert not residue_in_topology(new_topology, "HOH")


def test_remove_atoms_in_place(topology_with_ligand):
    new_topology = deepcopy(topology_with_ligand)
    new_topology.remove_atoms(
        [8, 9, 10], inplace=True
    )
    assert new_topology.n_atoms == 17
    assert new_topology.n_chains == 4
    assert not residue_in_topology(new_topology, "HOH")


def test_add_bonds_from_dict(topology_2_chains):
    topology_2_chains.add_bonds_from_dict(
        {0: [1, 2],
         1: [2, 3],
         }
    )
    assert topology_2_chains.n_bonds == 4


def test_non_hyd_indices(topology_with_hydrogen):
    indices = topology_with_hydrogen.non_hyd_indices()
    expected = [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]
    assert indices == expected


def test_get_atoms_residues(topology_2_chains):
    atoms = np.array([0, 1, 6])
    residues = topology_2_chains.get_atoms_residues(atoms)
    assert residues == [0, 2]


def test_get_atoms_residues_result_is_sorted(topology_2_chains):
    atoms = np.array([6, 1, 0])
    residues = topology_2_chains.get_atoms_residues(atoms)
    assert residues == [0, 2]


def test_atoms_residues_only_gets_aminoacids_and_ligands(topology_with_ligand):
    indices = [0, 4, 8, 11, 16]
    residues = topology_with_ligand.get_atoms_residues(indices)
    assert residues == [0, 1, 4]


def test_get_residues_atoms(topology_2_chains):
    residues = [1, 2]
    atoms = topology_2_chains.get_residues_atoms(residues)
    assert atoms == [3, 4, 5, 6, 7, 8]


def test_residues_subset(topology_2_chains):
    subset = topology_2_chains.residues_subset([1, 2])
    assert subset.n_chains == 2
    assert subset.n_atoms == 6
