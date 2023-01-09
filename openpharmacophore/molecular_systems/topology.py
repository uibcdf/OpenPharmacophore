import mdtraj
import mdtraj.core.element as mdt_element
import pyunitwizard as puw


class Topology:
    """ Wrapper for mdtraj's topology object.
    """
    def __init__(self, topology=None):
        if topology is None:
            self.top = mdtraj.Topology()
        else:
            self.top = topology

    @property
    def n_chains(self) -> int:
        return self.top.n_chains

    @property
    def n_residues(self) -> int:
        return self.top.n_residues

    @property
    def n_atoms(self) -> int:
        return self.top.n_atoms

    def set_num_chains(self, num):
        """ Set the number of chains in an empty topology.

            Parameters
            ----------
            num : int
                New number of chains
        """
        if self.n_chains == 0:
            for ii in range(num):
                self.top.add_chain()
        else:
            raise ValueError(f"Topology already has {self.n_chains} chains")

    def add_residue(self, res_name, chain):
        """ Add a new residue to a chain.

            Parameters
            ----------
            res_name : str
                Name of the residue

            chain : int
                Index of the chain to which the residue will be added

            Returns
            -------
            int
                The index of the new residue.
        """
        chain = self.top.chain(chain)
        self.top.add_residue(res_name, chain)
        return self.n_residues - 1

    def add_atom(self, name, symbol, residue):
        """ Add a new atom to the topology.

            Parameters
            ----------
            name : str
                Name of the atom

            symbol : str
                Symbol of the atom

            residue : int
                Index of the residue to which the atom will be added

            Returns
            -------
            int
                The index of the new atom
        """
        residue = self.top.residue(residue)
        element = mdt_element.get_by_symbol(symbol)
        self.top.add_atom(name, element, residue)
        return self.n_atoms - 1

    def add_atoms_to_chain(self, res_atoms, chain):
        """ Add multiple atoms to a chain.

            Parameters
            ----------
            res_atoms : dict[str, list[tuple[str, str]]
                Dictionary with the atoms names and symbols of each
                new residue.

            chain : int
                The index of the chain
        """
        for res_name, atoms in res_atoms.items():
            res_ind = self.add_residue(res_name, chain)
            for name, symbol in atoms:
                self.add_atom(name, symbol, res_ind)

    def has_hydrogens(self):
        """ Returns true if there is at least one hydrogen
            atom in the topology.

            Returns
            -------
            bool
        """
        for atom in self.top.atoms:
            if atom.element.symbol == "H":
                return True
        return False


def create_topology(traj_file, topology_file=None):
    """ Create a topology object.

        Parameters
        ----------
        traj_file : str
            Path to a trajectory file or pdb. If the trajectory file
            does not contain topology a topology file is needed.

        topology_file : str, optional
            File with the topology.

        Returns
        -------
        topology : Topology
            The topology

        coords : Quantity
            The coordinates of the system in angstroms.
    """
    if topology_file is not None:
        traj = mdtraj.load(traj_file, top=topology_file)
    else:
        traj = mdtraj.load(traj_file)

    topology = Topology(traj.topology)
    coords = puw.quantity(traj.xyz, "angstroms")
    return topology, coords
