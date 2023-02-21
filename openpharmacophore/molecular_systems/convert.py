import mdtraj as mdt
from mdtraj.utils.unit import in_units_of
from mdtraj.formats.pdb import PDBTrajectoryFile
from rdkit import Chem
import pyunitwizard as puw

import io

from openpharmacophore.molecular_systems.topology import Topology, CHAIN_NAMES
from openpharmacophore.constants import QuantityLike


class PDBStringIO(PDBTrajectoryFile):
    """ Class to write pdb files to a text stream or StringIO
    """
    def __init__(self):
        self._open = False
        self._file = None
        self._topology = None
        self._positions = None
        self._mode = "w"
        self._last_topology = None
        self._standard_names = True

        self._header_written = False
        self._footer_written = False
        self._file = io.StringIO()

        self._open = True

    @property
    def sio(self) -> io.StringIO:
        return self._file


def topology_to_mol(topology, coords, remove_hyd=True):
    """ Create a molecule from a topology object and an array of coordinates

           Parameters
           ----------
           topology : Topology
               Topology of the ligand.

           coords: np.ndarray
               Array of shape (n_atoms, 3) in nanometers.

           remove_hyd : bool
               Whether to remove the hydrogens from the molecule.

           Returns
           -------
           mol: rdkit.Chem.Mol

   """
    with PDBStringIO() as pdb:
        pdb.write(
            in_units_of(coords, mdt.Trajectory._distance_unit, pdb.distance_unit),
            topology.top,
            modelIndex=0,
            bfactors=None
        )

        mol = Chem.MolFromPDBBlock(pdb.sio.getvalue(), removeHs=remove_hyd)
    assert mol is not None, "Failed to create molecule"
    return mol


def _get_number_of_chains(mol):
    """ Count the number of chains in a rdkit molecule.

        Parameters
        ----------
        mol : rdkit.Chem.Mol

        Returns
        -------
        int
    """
    # Assumes that the chains are in order
    n_atoms = mol.GetNumAtoms()
    return CHAIN_NAMES.index(
        mol.GetAtomWithIdx(n_atoms - 1).GetPDBResidueInfo().GetChainId()
    ) + 1


def _add_residues(topology, mol, residue_map):
    """ Add residues to the topology.

        Parameters
        ----------
        topology : Topology

        mol : rdkit.Chem.Mol

        residue_map : dict[int, int]

    """
    res_ind = 0

    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        res_num = info.GetResidueNumber()

        if res_num not in residue_map:
            topology.add_residue(info.GetResidueName(),
                                 CHAIN_NAMES.index(info.GetChainId()))
            residue_map[res_num] = res_ind
            res_ind += 1


def _add_atoms(topology, mol, residue_map):
    """ Add atoms to the topology.

        Parameters
        ----------
        topology : Topology

        mol : rdkit.Chem.Mol

        residue_map : dict[int, int]

    """
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        topology.add_atom(
            name=info.GetName(),
            symbol=atom.GetSymbol(),
            residue=residue_map[info.GetResidueNumber()]
        )


def _add_bonds(topology, mol):
    """ Add bonds to the topology.

        Parameters
        ----------
        topology : Topology

        mol : rdkit.Chem.Mol

    """
    for bond in mol.GetBonds():
        at_1 = topology.get_atom(bond.GetBeginAtomIdx())
        at_2 = topology.get_atom(bond.GetEndAtomIdx())
        topology.add_bond(at_1, at_2)


def mol_to_topology(mol):
    """ Convert an rdkit molecule to a topology.

        Parameters
        ----------
        mol : rdkit.Chem.Mol

        Returns
        -------
        Topology

        Raises
        ------
        ValueError
            If the molecule does not contain pdb residue info.
    """
    topology = Topology()
    if mol.GetNumAtoms() == 0:
        return topology

    assert mol.GetAtomWithIdx(0).GetPDBResidueInfo() is not None, "Molecule does not contain PDB residue info"

    # residue numbers may not start from 1, so we must create a map
    residue_map = {}

    topology.set_num_chains(_get_number_of_chains(mol))
    _add_residues(topology, mol, residue_map)
    _add_atoms(topology, mol, residue_map)
    _add_bonds(topology, mol)
    return topology


def create_traj(coords, topology):
    """ Create a mdtraj Trajectory from an array of coordinates and
        a topology.

        Parameters
        ----------
        coords : QuantityLike

        topology : Topology

        Returns
        -------
        mdtraj.Trajectory
    """
    return mdt.Trajectory(
        xyz=puw.get_value(coords, "nanometers"),
        topology=topology.top
    )
