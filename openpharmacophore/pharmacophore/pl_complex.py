from .._private_tools import exceptions as exc
from ..data import pdb_to_smi
import numpy as np
import mdtraj as mdt
import rdkit.Chem.AllChem as Chem
import pyunitwizard as puw
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
            raise exc.NoLigandError

        self.coords = self.traj.xyz
        self.mol_graph = Chem.MolFromPDBFile(file_path)

        self._receptor_indices = []

        self._lig_indices = []
        self._ligand = None
        self._ligand_id = ""

    @property
    def ligand_ids(self):
        return self._ligand_ids

    @property
    def ligand(self):
        return self._ligand

    def set_ligand(self, lig_id):
        """ Set the ligand that will be used to do all computations.

            Parameters
            ----------
            lig_id : str
                Ligand id.
        """
        if lig_id not in self._ligand_ids:
            raise ValueError("Invalid ligand id")
        self._ligand_id = lig_id

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

    def _ligand_and_receptor_indices(self):
        """ Get the indices of the ligand with the given id and
            those of the receptor.

            Parameters
            ----------
            lig_id : str
                Ligand id.
        """
        ligand, chain = self._ligand_id.split(":")
        chain_index = self.chain_names.index(chain)
        for atom in self.topology.atoms:
            if atom.residue.name == ligand and atom.residue.chain.index == chain_index:
                self._lig_indices.append(atom.index)
            else:
                self._receptor_indices.append(atom.index)

    def _ligand_to_mol(self):
        """ Extract the ligand from the trajectory and create and rdkit mol.
        """
        if len(self._lig_indices) == 0:
            raise exc.NoLigandIndicesError

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

    def _remove_ligand(self):
        """ Remove a ligand from the trajectory.
        """
        if len(self._lig_indices) == 0:
            raise exc.NoLigandIndicesError

        self.traj = self.traj.atom_slice(self._receptor_indices)
        self.topology = self.traj.topology

    def has_hydrogens(self):
        """ Returns true if the topology contains hydrogens.
        """
        for atom in self.topology.atoms:
            if atom.element.symbol == "H":
                return True
        return False

    def fix_ligand(self, smiles=""):
        """ Add hydrogens to the ligand and correct its bond orders.

            If a smiles is not given, a smiles for the ligand will
            be searched for. In case it is not found this will raise
            an error.

            Parameters
            ----------
            smiles : str, optional
                The smiles of the ligand.

            Raises
            ------
            SmilesNotFoundError
                If a smiles for the ligand is not found.
        """
        if not smiles:
            smiles = self._pdb_id_to_smi(self._ligand_id)

        template = Chem.MolFromSmiles(smiles)
        if self._ligand.GetNumAtoms() != template.GetNumAtoms():
            raise exc.DifferentNumAtomsError(
                self._ligand.GetNumAtoms(), template.GetNumAtoms()
            )

        self._ligand = Chem.AssignBondOrdersFromTemplate(template, self._ligand)

        # TODO: Add hydrogens

    @staticmethod
    def _pdb_id_to_smi(pdb_id):
        """ Get the smiles of a ligand from its pdb id.

            Parameters
            ----------
            pdb_id: str
                The pdb id of the ligand

            Returns
            -------
            str
                The smiles of the ligand.
        """
        lig_name = pdb_id.split(":")[0]
        if lig_name == "UNL":
            raise exc.SmilesNotFoundError(lig_name)

        with open(pdb_to_smi) as fp:
            lines = fp.readlines()

        pdb_id_mapper = {}
        for line in lines:
            lig_id, smi = line.split()
            pdb_id_mapper[lig_id] = smi

        try:
            return pdb_id_mapper[lig_name]
        except KeyError:
            raise exc.SmilesNotFoundError(lig_name)

    @staticmethod
    def _modeller_to_trajectory(modeller):
        """ Convert an openmm.Modeller to a mdtraj.Trajectory.

            Parameters
            ----------
            modeller : openmm.Modeller

            Returns
            -------
            mdtraj.Trajectory
        """
        positions = modeller.getPositions()
        n_atoms = len(positions)
        coords = np.zeros((1, n_atoms, 3))

        for jj in range(n_atoms):
            coords[0, jj, :] = puw.get_value(
                positions[jj], to_unit="nanometers")

        topology = mdt.Topology.from_openmm(modeller.getTopology())
        return mdt.Trajectory(coords, topology)
