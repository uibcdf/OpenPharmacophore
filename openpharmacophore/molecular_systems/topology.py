import mdtraj
import mdtraj.core.element as mdt_element
import pyunitwizard as puw

SOLVENT_AND_IONS = frozenset(
    ['118', '119', '1AL', '1CU', '2FK', '2HP', '2OF',
     '3CO', '3MT', '3NI', '3OF', '4MO', '543', '6MO', 'ACT', 'AG', 'AL', 'ALF',
     'ATH', 'AU', 'AU3', 'AUC', 'AZI', 'Ag', 'BA', 'BAR', 'BCT', 'BEF', 'BF4',
     'BO4', 'BR', 'BS3', 'BSY', 'Be', 'CA', 'CA+2', 'Ca+2', 'CAC', 'CAD', 'CAL',
     'CD', 'CD1', 'CD3', 'CD5', 'CE', 'CES', 'CHT', 'CL', 'CL-', 'CLA', 'Cl-', 'CO',
     'CO3', 'CO5', 'CON', 'CR', 'CS', 'CSB', 'CU', 'CU1', 'CU3', 'CUA', 'CUZ',
     'CYN', 'Cl-', 'Cr', 'DME', 'DMI', 'DSC', 'DTI', 'DY', 'E4N', 'EDR', 'EMC',
     'ER3', 'EU', 'EU3', 'F', 'FE', 'FE2', 'FPO', 'GA', 'GD3', 'GEP', 'HAI', 'HG',
     'HGC', 'HOH', 'IN', 'IOD', 'ION', 'IR', 'IR3', 'IRI', 'IUM', 'K', 'K+', 'KO4',
     'LA', 'LCO', 'LCP', 'LI', 'LIT', 'LU', 'MAC', 'MG', 'MH2', 'MH3', 'MLI', 'MMC',
     'MN', 'MN3', 'MN5', 'MN6', 'MO1', 'MO2', 'MO3', 'MO4', 'MO5', 'MO6', 'MOO',
     'MOS', 'MOW', 'MW1', 'MW2', 'MW3', 'NA', 'NA+2', 'NA2', 'NA5', 'NA6', 'NAO',
     'NAW', 'Na+2', 'NET', 'NH4', 'NI', 'NI1', 'NI2', 'NI3', 'NO2', 'NO3', 'NRU',
     'Na+', 'O4M', 'OAA', 'OC1', 'OC2', 'OC3', 'OC4', 'OC5', 'OC6', 'OC7', 'OC8',
     'OCL', 'OCM', 'OCN', 'OCO', 'OF1', 'OF2', 'OF3', 'OH', 'OS', 'OS4', 'OXL',
     'PB', 'PBM', 'PD', 'PER', 'PI', 'PO3', 'PO4', 'POT', 'PR', 'PT', 'PT4', 'PTN',
     'RB', 'RH3', 'RHD', 'RU', 'RUB', 'Ra', 'SB', 'SCN', 'SE4', 'SEK', 'SM', 'SMO',
     'SO3', 'SO4', 'SOD', 'SR', 'Sm', 'Sn', 'T1A', 'TB', 'TBA', 'TCN', 'TEA', 'THE',
     'TL', 'TMA', 'TRA', 'UNX', 'V', 'V2+', 'VN3', 'VO4', 'W', 'WO5', 'Y1', 'YB',
     'YB2', 'YH', 'YT3', 'ZN', 'ZN2', 'ZN3', 'ZNA', 'ZNO', 'ZO3']
)

CHAIN_NAMES = [
    "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K",
    "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V",
    "W", "X", "Y", "Z",
]


class Topology:
    """ Wrapper for mdtraj's Topology object.
    """

    def __init__(self, topology=None):
        if topology is None:
            self.top = mdtraj.Topology()
        else:
            self.top = topology  # type: mdtraj.Topology

    @property
    def n_chains(self) -> int:
        return self.top.n_chains

    @property
    def n_residues(self) -> int:
        return self.top.n_residues

    @property
    def n_atoms(self) -> int:
        return self.top.n_atoms

    @property
    def n_bonds(self) -> int:
        return self.top.n_bonds

    @classmethod
    def from_openmm(cls, openmm_top):
        """ Create a topology from an openmm.Topology.
        """
        return cls(topology=mdtraj.Topology.from_openmm(openmm_top))

    def iter_atoms(self):
        """ Iterator of the atoms of the topology
        """
        return self.top.atoms

    def set_num_chains(self, num):
        """ Set the number of chains in an empty topology.

            Parameters
            ----------
            num : int
                New number of chains
        """
        if self.n_chains == 0:
            for ii in range(num):
                self.add_chain()
        else:
            raise ValueError(f"Topology already has {self.n_chains} chains")

    def add_chain(self):
        """ Adda new chain.

            Returns
            -------
            int
                Index of the new chain
        """
        self.top.add_chain()
        return self.n_chains - 1

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

    def add_bond(self, at_ind_1, at_ind_2):
        """ Add a new atom to the topology

            Parameters
            ----------
            at_ind_1 : int
                Index of the first atom.

            at_ind_2
                Index of the second atom.
        """
        self.top.add_bond(at_ind_1, at_ind_2)

    def get_atom(self, index):
        """ Get the atom with the given index.
        """
        return self.top.atom(index)

    def get_residue(self, index):
        """ Get the residue with given index.
        """
        return self.top.residue(index)

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

    @staticmethod
    def _is_ligand_atom(atom):
        """ Check if an atom belongs to a ligand.

            Parameters
            ----------
            atom : mdtraj.Atom

            Returns
            -------
            bool
        """
        return not atom.residue.is_water and not atom.residue.is_protein \
            and atom.residue.name not in SOLVENT_AND_IONS

    def has_ligands(self):
        """ Returns true if there are any ligands in the topology.

            Returns
            -------
            bool
        """
        for atom in self.top.atoms:
            if self._is_ligand_atom(atom):
                return True
        return False

    def ligand_ids(self):
        """ Returns the ligand ids of the ligands in the protein if there
            are any.

            Returns
            -------
            ligands : list[str]
                List with ligand ids. An empty list if there are no ligands
        """
        ligands = set()
        for atom in self.top.atoms:
            if self._is_ligand_atom(atom):
                chain = atom.residue.chain.index
                ligand_id = atom.residue.name + ":" + CHAIN_NAMES[chain]
                ligands.add(ligand_id)

        return list(ligands)

    def get_residue_indices(self, res_name, chain):
        """ Get the indices of hte residue with given name and chain

            Parameters
            ----------
            res_name : str
                Name of the residue

            chain : str
                Name of the chain

            Returns
            -------
            indices : list[int]

        """
        chain_index = CHAIN_NAMES.index(chain)
        indices = [
            a.index for a in self.top.atoms
            if a.residue.name == res_name and a.residue.chain.index == chain_index
        ]
        return indices

    def remove_atoms(self, indices, inplace=False):
        """ Remove atoms from the topology.

            Parameters
            ----------
            indices : list[int]
                Indices of the atoms that will be removed

            inplace : bool, default=False
                Whether to perform the operation in place.

            Returns
            -------
            Topology
                A new topology with the specified atoms removed, if inplace
                its false. Otherwise, it does not return anything.

        """
        # indices of atoms that will not be removed
        atoms = [
            ii for ii in range(0, self.n_atoms) if ii not in indices
        ]
        if inplace:
            self.top = self.top.subset(atoms)
        else:
            return Topology(self.top.subset(atoms))

    def add_bonds_from_dict(self, bonds_dict):
        """ Add bonds to the topology

            Parameters
            ----------
            bonds_dict: dict [int, list[int]]
                Each entry represents an atom index (key) to which other atoms
                are bonded (values).

        """
        for atom_ind, neighbors_ind in bonds_dict.items():
            atom = self.top.atom(atom_ind)
            for nbr in neighbors_ind:
                self.top.add_bond(atom, self.top.atom(nbr))

    def subset(self, atom_indices):
        """ Returns a subset of the topology.

            Parameters
            ----------
            atom_indices : list[int]

            Returns
            -------
            Topology
        """
        return Topology(self.top.subset(atom_indices))

    def residues_subset(self, residues):
        """ Returns a subset of the topology.

            Parameters
            ----------
            residues : list[int]

            Returns
            -------
            Topology

        """
        return self.subset(self.get_residues_atoms(residues))

    def non_hyd_indices(self):
        """ Returns the indices of the atoms that are not hydrogen.

            Returns
            -------
            list[int]
        """
        return [
            a.index for a in self.top.atoms if a.element.symbol != "H"
        ]

    def to_openmm(self):
        """ Convert to openmm topology. """
        return self.top.to_openmm()

    def get_bs_residues(self, atoms):
        """ Get the indices of the residues in the binding site.

            Parameters
            ----------
            atoms : ArrayLike[int]
                An array like object with the atoms indices

            Returns
            -------
            list[int]
                The residues indices in ascending order

            int
                Index of the ligand. None if there is no ligand
        """
        residues = set()
        ligand = None
        for at_ind in atoms:
            atom = self.top.atom(at_ind)
            if atom.residue.is_protein:
                residues.add(atom.residue.index)
            elif self._is_ligand_atom(atom):
                ligand = atom.residue.index

        return sorted(residues), ligand

    def get_residues_atoms(self, residues):
        """ Get the indices of the atoms which comprise the given residues.

            Parameters
            ----------
            residues : list[int]
                An array like object with the residues indices. It must be sorted

            Returns
            -------
            atoms : list[int]
                A list with atoms indices
        """
        atoms = []
        for res in residues:
            for at in self.top.residue(res).atoms:
                atoms.append(at.index)
        return atoms

    def join(self, other):
        """ Join two topologies.

            Parameters
            ----------
            other : Topology

            Returns
            -------
            Topology
        """
        return Topology(self.top.join(other.top))

    def residue_has_hyd(self, residue_ind):
        """ Check if a residue has hydrogens

        Parameters
        ----------
        residue_ind : int

        Returns
        -------
        bool

        """
        for atom in self.get_residue(residue_ind).atoms:
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
    coords = puw.quantity(traj.xyz, "nanometers")
    return topology, coords
