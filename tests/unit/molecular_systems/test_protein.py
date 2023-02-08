import pyunitwizard as puw
import numpy as np
from copy import deepcopy

from openpharmacophore import Protein


def test_remove_ligand(topology_with_ligand):
    topology = deepcopy(topology_with_ligand)
    coords = puw.quantity(
        np.ones((2, topology.n_atoms, 3)),
        "nanometers"
    )
    protein = Protein(topology, coords)
    n_atoms = protein.n_atoms

    protein.remove_ligand("EST:B")
    assert protein.n_atoms == n_atoms - 4
    assert protein.coords.shape == (2, n_atoms - 4, 3)

    expected_coords = np.ones((2, n_atoms - 4, 3))
    assert np.all(puw.get_value(protein.coords, "nanometers")
                  == expected_coords)


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
