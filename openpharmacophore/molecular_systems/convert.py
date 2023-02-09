import mdtraj as mdt
from mdtraj.utils.unit import in_units_of
from mdtraj.formats.pdb import PDBTrajectoryFile
from rdkit import Chem
import pyunitwizard as puw

import io

from openpharmacophore.molecular_systems import Topology
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

    # Add first chain and residue
    atom = mol.GetAtomWithIdx(0)
    info = atom.GetPDBResidueInfo()
    if info is None:
        # Assume that if the first atom has residue info all of them have
        raise ValueError("Molecule does not contain PDB residue info")
    prev_chain = info.GetChainId()
    prev_res = info.GetResidueName()

    chain = topology.add_chain()
    residue = topology.add_residue(info.GetResidueName(), chain)

    # Add the rest of chains and residues as well as atoms
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        chain_id = info.GetChainId()
        res = info.GetResidueName()
        if chain_id != prev_chain:
            chain = topology.add_chain()
            prev_chain = chain_id
        if prev_res != res:
            residue = topology.add_residue(res, chain)
            prev_res = res

        topology.add_atom(info.GetName(), atom.GetSymbol(), residue)

    for bond in mol.GetBonds():
        at_1 = topology.get_atom(bond.GetBeginAtomIdx())
        at_2 = topology.get_atom(bond.GetEndAtomIdx())
        topology.add_bond(at_1, at_2)

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
