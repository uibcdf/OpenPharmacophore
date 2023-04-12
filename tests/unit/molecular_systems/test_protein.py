import pyunitwizard as puw
import numpy as np
from copy import deepcopy

from openpharmacophore.molecular_systems import Protein, Topology


def protein_with_ligand(topology_with_ligand):
    topology = topology_with_ligand
    coords = puw.quantity(
        np.ones((2, topology.n_atoms, 3)),
        "nanometers"
    )
    return Protein(topology, coords)


def test_remove_ligand(topology_with_ligand):
    protein = protein_with_ligand(topology_with_ligand)
    n_atoms = protein.n_atoms

    protein.remove_ligand("EST:B")
    assert protein.n_atoms == n_atoms - 4
    assert protein.coords.shape == (2, n_atoms - 4, 3)

    expected_coords = np.ones((2, n_atoms - 4, 3))
    assert np.all(
        puw.get_value(protein.coords, "nanometers") == expected_coords
    )


def test_atoms_at_distance(protein_4_residues):
    centroid = puw.quantity(np.zeros((1, 3)), "nanometers")
    max_dist = puw.quantity(3.0, "nanometers")

    indices = protein_4_residues.atoms_at_distance(0, centroid, max_dist)
    expected = np.array([1, 3, 8])
    assert np.all(indices == expected)


def test_slice_protein_all_frames(protein_4_residues):
    atoms = [3, 4, 5, 6, 7, 8]
    sliced_protein = protein_4_residues.slice(atoms)

    assert sliced_protein.n_atoms == 6
    assert sliced_protein.n_residues == 2
    assert sliced_protein.coords.shape == (1, 6, 3)

    expected = puw.quantity(np.array([[
        [1., 1., 1.],
        [4., 4., 4.],
        [4., 4., 4.],
        [4., 4., 4.],
        [4., 4., 4.],
        [1., 1., 1.],
    ]]), "nanometers")

    assert np.all(expected == sliced_protein.coords)


def test_slice_protein_single_frame(protein_4_residues):
    atoms = [3, 4, 5, 6, 7, 8]
    sliced_protein = protein_4_residues.slice(atoms, frame=0)

    assert sliced_protein.n_atoms == 6
    assert sliced_protein.n_residues == 2
    assert sliced_protein.coords.shape == (1, 6, 3)

    expected = puw.quantity(np.array([[
        [1., 1., 1.],
        [4., 4., 4.],
        [4., 4., 4.],
        [4., 4., 4.],
        [4., 4., 4.],
        [1., 1., 1.],
    ]]), "nanometers")

    assert np.all(expected == sliced_protein.coords)


def test_concatenate_to_protein(protein_4_residues, estradiol_topology, estradiol_coords):
    protein = deepcopy(protein_4_residues)
    start_atoms = protein.n_atoms
    estradiol_atoms = estradiol_topology.n_atoms

    protein.concatenate(estradiol_topology, estradiol_coords)
    assert protein.n_atoms == start_atoms + estradiol_atoms
    assert protein.coords.shape == (1, start_atoms + estradiol_atoms, 3)


def glycine():
    topology = Topology()
    topology.add_chain()
    topology.add_atoms_to_chain({
        "GLY": [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("O", "O")]
    }, 0)
    coords = puw.quantity(np.array([[
        [1.9310, 0.0900, -0.0340],
        [0.7610, -0.7990, -0.0080],
        [-0.4980, 0.0290, -0.0050],
        [-0.4290, 1.2350, -0.0230],
        [-1.6970, -0.5740, 0.0180],
    ]]), "angstroms")
    return Protein(topology, coords)


def test_add_hydrogens_to_protein():
    protein = glycine()
    protein.add_hydrogens()

    assert protein.has_hydrogens
    assert protein.n_atoms == 10


def test_extract_chain(topology_2_chains):
    coords = np.arange(1, topology_2_chains.n_atoms * 3 + 1).reshape((1, -1, 3))
    coords = puw.quantity(coords, "angstroms")

    prot = Protein(topology_2_chains, coords)
    prot.extract_chain("A")

    assert prot.n_atoms == 6
    assert prot.n_chains == 1

    expected_coords = coords[:, 0:6, :]
    assert np.all(prot.coords == expected_coords)


def test_remove_all_ligands(topology_with_ligand):
    prot = protein_with_ligand(topology_with_ligand)
    prot.remove_all_ligands()

    assert not prot.has_ligands()
    assert len(prot.ligand_ids()) == 0
    assert prot.n_atoms == 12
    assert prot.n_chains == 3


def test_remove_solvent_and_ions(topology_with_ligand):
    prot = protein_with_ligand(topology_with_ligand)
    prot.remove_solvent_and_ions()

    assert prot.n_atoms == 12
    assert prot.n_chains == 3
    assert not prot.has_solvent_or_ions()
