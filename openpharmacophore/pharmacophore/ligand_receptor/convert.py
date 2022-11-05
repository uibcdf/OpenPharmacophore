import mdtraj as mdt


def mol_to_topology(mol):
    """ Convert a rdkit molecule to a mdtraj topology.

        Parameters
        ----------
        mol : rdkit.Chem.Mol

        Returns
        -------
        mdtraj.Topology
    """
    pass


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
    pass
