import mdtraj as mdt
import mdtraj.core.element as mdt_element
import numpy as np


chain_names = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K",
               "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V",
               "W", "X", "Y", "Z"]


class NoConformersError(ValueError):
    """ Exception raised when a molecule contains no conformers.
    """
    pass


def mol_to_topology(mol):
    """ Convert a rdkit molecule to a mdtraj topology.

        Parameters
        ----------
        mol : rdkit.Chem.Mol

        Returns
        -------
        mdtraj.Topology

        Raises
        ------
        ValueError
            If the molecule does not contain pdb residue info.
    """
    topology = mdt.Topology()
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
        element = mdt_element.get_by_symbol(atom.GetSymbol())
        topology.add_atom(info.GetName(), element, residue)

    for bond in mol.GetBonds():
        at_1 = topology.atom(bond.GetBeginAtomIdx())
        at_2 = topology.atom(bond.GetEndAtomIdx())
        topology.add_bond(at_1, at_2)

    return topology


def mol_to_traj(mol):
    """ Convert a rdkit molecule to a mdtraj trajectory.

        A frame will be created for each conformer.

        Parameters
        ----------
        mol : rdkit.Chem.Mol

        Returns
        -------
        mdtraj.Trajectory

        Raises
        ------
        NoConformersError
            If the molecule has no conformers.
        """
    if mol.GetNumConformers() == 0:
        raise NoConformersError

    topology = mol_to_topology(mol)
    n_confs = mol.GetNumConformers()
    n_atoms = mol.GetNumAtoms()
    coords = np.zeros((n_confs, n_atoms, 3))

    for ii, conf in enumerate(mol.GetConformers()):
        for jj in range(n_atoms):
            pos = conf.GetAtomPosition(jj)
            # Convert to nanometers
            coords[ii, jj, 0] = pos.x / 10
            coords[ii, jj, 1] = pos.y / 10
            coords[ii, jj, 2] = pos.z / 10

    return mdt.Trajectory(coords, topology)
