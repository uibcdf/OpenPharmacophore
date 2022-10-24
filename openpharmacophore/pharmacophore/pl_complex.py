from .._private_tools.exceptions import NoLigandError, NoLigandIndicesError
import mdtraj as mdt
import rdkit.Chem.AllChem as Chem
import tempfile


class PLComplex:
    """ Class to store protein-ligand complexes and compute their interactions.

        Parameters
        ----------
        file_path: str
            Path to a pdb file.

    """
    chain_names = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K",
                   "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V",
                   "W", "X", "Y", "Z"]

    def __init__(self, file_path):
        self.traj = mdt.load(file_path)
        self.topology = self.traj.topology

        self._ligand_ids = self.find_ligands()
        if len(self._ligand_ids) == 0:
            raise NoLigandError

        self.coords = self.traj.xyz
        self.mol_graph = Chem.MolFromPDBFile(file_path)

        self._lig_indices = []
        self._ligand = None

    @property
    def ligand_ids(self):
        return self._ligand_ids

    @property
    def ligand(self):
        return self._ligand

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
        if not atom.residue.is_water and not atom.residue.is_protein \
                and atom.residue.n_atoms > 4:  # Discard ions
            return True
        return False

    def find_ligands(self):
        """ Returns the ligand ids of the ligands in the complex.

            Returns
            -------
            ligands : list[str]
                List with ligand ids.
        """
        ligands = []
        for atom in self.topology.atoms:
            if self._is_ligand_atom(atom):
                chain = atom.residue.chain.index
                ligand_id = atom.residue.name + ":" + self.chain_names[chain]
                if ligand_id not in ligands:
                    ligands.append(ligand_id)

        return ligands

    def _ligand_atom_indices(self, lig_id):
        """ Get the indices of the ligand with the given id.

            Parameters
            ----------
            lig_id : str
                Ligand id.
        """
        ligand, chain = lig_id.split(":")
        chain_index = self.chain_names.index(chain)
        for atom in self.topology.atoms:
            if atom.residue.name == ligand and atom.residue.chain.index == chain_index:
                self._lig_indices.append(atom.index)

    def _ligand_to_mol(self):
        """ Extract the ligand from the trajectory and create and rdkit mol.
        """
        if len(self._lig_indices) == 0:
            raise NoLigandIndicesError

        # TODO: implement the conversion from trajectory to rdkit mol without using
        #  a file.
        lig_traj = self.traj.atom_slice(self._lig_indices)
        pdb_file = tempfile.NamedTemporaryFile()
        lig_traj.save_pdb(pdb_file.name)

        pdb_file.seek(0)
        mol = Chem.MolFromPDBFile(pdb_file.name)
        assert mol is not None
        pdb_file.close()

        self._ligand = mol
