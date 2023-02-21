from openpharmacophore.molecular_systems import Topology
import openpharmacophore.molecular_systems.convert as convert
from rdkit.Chem import AllChem as Chem
import pytest
import numpy as np


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
        topology : Topology
    """
    assert mol.GetNumAtoms() == topology.n_atoms, "Num atoms"
    assert mol.GetNumBonds() == topology.n_bonds, "Num bonds"

    elements_mol = sorted([a.GetSymbol() for a in mol.GetAtoms()])
    elements_top = sorted([a.element.symbol for a in topology.iter_atoms()])
    assert elements_mol == elements_top, "Element symbols"


def test_mol_with_no_residue_info_raises_error():
    mol = Chem.MolFromSmiles("CC(C(C(=O)O)N)O")
    with pytest.raises(AssertionError):
        convert.mol_to_topology(mol)


def test_mol_to_topology():
    estradiol = Chem.MolFromPDBBlock("""
MODEL        0
ATOM   5944  C1  EST A 600     104.106  17.203  24.775  1.00  0.00           C  
ATOM   5945  C2  EST A 600     102.995  17.834  25.370  1.00  0.00           C  
ATOM   5946  C3  EST A 600     101.695  17.355  25.120  1.00  0.00           C  
ATOM   5947  O3  EST A 600     100.598  17.990  25.704  1.00  0.00           O  
ATOM   5948  C4  EST A 600     101.506  16.240  24.274  1.00  0.00           C  
ATOM   5949  C5  EST A 600     102.621  15.588  23.660  1.00  0.00           C  
ATOM   5950  C6  EST A 600     102.371  14.379  22.735  1.00  0.00           C  
ATOM   5951  C7  EST A 600     103.644  13.753  22.086  1.00  0.00           C  
ATOM   5952  C8  EST A 600     104.898  13.873  22.953  1.00  0.00           C  
ATOM   5953  C9  EST A 600     105.178  15.388  23.261  1.00  0.00           C  
ATOM   5954  C10 EST A 600     103.957  16.078  23.918  1.00  0.00           C  
ATOM   5955  C11 EST A 600     106.462  15.459  24.125  1.00  0.00           C  
ATOM   5956  C12 EST A 600     107.711  14.803  23.508  1.00  0.00           C  
ATOM   5957  C13 EST A 600     107.463  13.343  23.124  1.00  0.00           C  
ATOM   5958  C14 EST A 600     106.170  13.270  22.242  1.00  0.00           C  
ATOM   5959  C15 EST A 600     106.228  11.821  21.792  1.00  0.00           C  
ATOM   5960  C16 EST A 600     107.701  11.713  21.263  1.00  0.00           C  
ATOM   5961  C17 EST A 600     108.494  12.719  22.135  1.00  0.00           C  
ATOM   5962  O17 EST A 600     109.610  12.027  22.746  1.00  0.00           O  
ATOM   5963  C18 EST A 600     107.379  12.449  24.419  1.00  0.00           C  
TER    5964      EST A 600
ENDMDL
CONECT    1    2   11
CONECT    2    1    3
CONECT    3    2    4    5
CONECT    4    3
CONECT    5    3    6
CONECT    6    5    7   11
CONECT    7    6    8
CONECT    8    7    9
CONECT    9    8   10   15
CONECT   10    9   11   12
CONECT   11    1    6   10
CONECT   12   10   13
CONECT   13   12   14
CONECT   14   13   15   18   20
CONECT   15    9   14   16
CONECT   16   15   17
CONECT   17   16   18
CONECT   18   14   17   19
CONECT   19   18
CONECT   20   14
END
""")
    assert estradiol.GetNumAtoms() == 20
    est_topology = convert.mol_to_topology(estradiol)

    assert_mol_and_topology_equal(estradiol, est_topology)
    assert est_topology.get_residue(0).name == "EST"


def test_get_number_of_chains(dialanine):
    assert convert._get_number_of_chains(dialanine) == 1


def test_add_residues_when_converting_mol_to_topology(dialanine):
    topology = Topology()
    topology.set_num_chains(1)
    convert._add_residues(topology, dialanine, dict())

    assert topology.n_residues == 2


def test_mol_to_topology_does_not_create_additional_residues(dialanine):
    topology = convert.mol_to_topology(dialanine)

    assert_mol_and_topology_equal(dialanine, topology)
    assert topology.n_residues == 2


@pytest.fixture
def dipeptide_topology():
    topology = Topology()
    topology.set_num_chains(1)
    topology.add_atoms_to_chain({
        "THR": [
            ("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"),
            ("CB", "C"), ("OG1", "O"), ("CG2", "C"),
        ],
        "TYR": [
            ("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"),
            ("CB", "C"), ("CG", "C"), ("CD1", "C"), ("CD2", "C"),
            ("CE1", "C"), ("CE2", "C"), ("CZ", "C"), ("OH", "O"),
        ],
    }, 0)
    return topology


@pytest.fixture
def dipeptide_coords():
    return np.array([[
        [4.4235, 8.0308, 1.8419],
        [4.3549, 7.9243, 1.7706],
        [4.4528, 7.8252, 1.7077],
        [4.5699, 7.8559, 1.6853],
        [4.2608, 7.9792, 1.6611],
        [4.3375, 8.0468, 1.5608],
        [4.1586, 8.0762, 1.7208],
        [4.4030, 7.7052, 1.6814],
        [4.4799, 7.6021, 1.6156],
        [4.4189, 7.5782, 1.4791],
        [4.3007, 7.6033, 1.4583],
        [4.4721, 7.4725, 1.6954],
        [4.5253, 7.4836, 1.8355],
        [4.6598, 7.4629, 1.8621],
        [4.4412, 7.5147, 1.9416],
        [4.7098, 7.4729, 1.9906],
        [4.4895, 7.5245, 2.0707],
        [4.6246, 7.5036, 2.0946],
        [4.6748, 7.5126, 2.2224],
    ]])


def test_topology_to_mol(dipeptide_topology, dipeptide_coords):
    mol = convert.topology_to_mol(
        dipeptide_topology, dipeptide_coords[0], remove_hyd=True
    )
    assert mol.GetNumAtoms() == 19
