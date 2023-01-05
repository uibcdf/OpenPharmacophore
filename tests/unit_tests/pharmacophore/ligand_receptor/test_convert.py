import openpharmacophore.pharmacophore.ligand_receptor.convert as convert
from tst_data import threonine, estradiol
import numpy as np
import pytest
from rdkit.Chem import AllChem as Chem
from copy import deepcopy


def sort_bonds(bonds):
    sorted_bonds = []
    for bnd in bonds:
        if bnd[0] < bnd[1]:
            sorted_bonds.append((bnd[0], bnd[1]))
        else:
            sorted_bonds.append((bnd[1], bnd[0]))
    sorted_bonds.sort()
    return sorted_bonds


def assert_mol_and_topology_equal(mol, topology):
    """ Assert that a molecule and a topology are the same.

        Parameters
        ----------
        mol : rdkit.Chem.Mol
        topology : mdtraj.Topology
    """
    assert mol.GetNumAtoms() == topology.n_atoms, "Num atoms"
    assert mol.GetNumBonds() == topology.n_bonds, "Num bonds"

    elements_mol = [a.GetSymbol() for a in mol.GetAtoms()]
    elements_top = [a.element.symbol for a in topology.atoms]
    assert elements_mol == elements_top, "Element symbols"

    chains_mol = set()
    residues_mol = []
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        chains_mol.add(info.GetChainId())
        residues_mol.append(info.GetResidueName())
    assert len(chains_mol) == topology.n_chains, "Num chain"
    assert len(set(residues_mol)) == topology.n_residues, "Num residues"

    residues_top = [a.residue.name for a in topology.atoms]
    assert residues_mol == residues_top, "Different residues"

    bonds_mol = [
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
    ]
    bonds_mol = sort_bonds(bonds_mol)

    bonds_top = [
        (b[0].index, b[1].index) for b in topology.bonds
    ]
    bonds_top = sort_bonds(bonds_top)

    assert bonds_top == bonds_mol, "Different bonds"


def test_mol_with_no_residue_info_raises_error():
    mol = Chem.MolFromSmiles("CC(C(C(=O)O)N)O")
    with pytest.raises(ValueError):
        convert.mol_to_topology(mol)


def test_mol_with_single_residue_to_topology():
    thr_topology = convert.mol_to_topology(threonine)
    assert_mol_and_topology_equal(threonine, thr_topology)
    assert thr_topology.residue(0).name == "THR"

    est_topology = convert.mol_to_topology(estradiol)
    assert_mol_and_topology_equal(estradiol, est_topology)
    assert est_topology.residue(0).name == "EST"


def test_mol_with_no_conformers_raises_error():
    mol = deepcopy(threonine)
    mol.RemoveAllConformers()
    assert mol.GetNumConformers() == 0
    with pytest.raises(convert.NoConformersError):
        convert.mol_to_traj(mol)


def test_mol_to_traj_mol_with_single_conformer():
    thr_traj = convert.mol_to_traj(threonine)
    n_atoms = threonine.GetNumAtoms()
    expected_coords = np.array([
        [44.235,  80.308,  18.419],
        [43.549,  79.243,  17.706],
        [44.528,  78.252,  17.077],
        [45.699,  78.559,  16.853],
        [42.608,  79.792,  16.611],
        [43.375,  80.468,  15.608],
        [41.586,  80.762,  17.208],
    ]) / 10  # Mdtraj expects nanometers

    assert thr_traj.n_atoms == n_atoms
    assert thr_traj.xyz.shape == (1, n_atoms, 3)
    assert np.allclose(thr_traj.xyz, expected_coords)

    est_traj = convert.mol_to_traj(estradiol)
    n_atoms = estradiol.GetNumAtoms()
    assert est_traj.n_atoms == n_atoms
    assert est_traj.xyz.shape == (1, n_atoms, 3)


@pytest.fixture()
def est_hyd():
    """ Returns estradiol molecule with hydrogens. """
    template = Chem.MolFromSmiles("C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O")
    mol = Chem.AssignBondOrdersFromTemplate(template, estradiol)
    mol = Chem.AddHs(mol, addCoords=True, addResidueInfo=True)
    assert mol.GetNumAtoms() == 44
    return mol


def test_molecule_with_hydrogens_to_topology(est_hyd):
    est_topology = convert.mol_to_topology(est_hyd)
    assert_mol_and_topology_equal(est_hyd, est_topology)
