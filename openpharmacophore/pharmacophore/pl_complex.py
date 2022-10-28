from .._private_tools import exceptions as exc
from ..data import pdb_to_smi
from ..utils import maths
import mdtraj as mdt
from nglview import show_mdtraj
import numpy as np
from openmm.app import Modeller
import pyunitwizard as puw
import rdkit.Chem.AllChem as Chem
import tempfile
from matplotlib.colors import to_rgb


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

    BS_DIST_MAX = puw.quantity(0.85, "nanometers")

    smarts_ligand = {
        "aromatic ring": [
            "a1aaaa1",
            "a1aaaaa1"
        ],
        "hydrophobicity": [
            '[$([S]~[#6])&!$(S~[!#6])]',
            '[C&r3]1~[C&r3]~[C&r3]1',
            '[C&r4]1~[C&r4]~[C&r4]~[C&r4]1',
            '[C&r5]1~[C&r5]~[C&r5]~[C&r5]~[C&r5]1',
            '[C&r6]1~[C&r6]~[C&r6]~[C&r6]~[C&r6]~[C&r6]1',
            '[C&r7]1~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]1',
            '[C&r8]1~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]1',
            '[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]',
            '*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]',
            '[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,'
            'I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,'
            'CH2X3,CH1X2,F,Cl,Br,I]',
            '[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]',
            '[$([CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,'
            'CH0X2])])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2] ',
            '[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]',
        ],
        "negative charge": [
            '[$([-,-2,-3])&!$(*[+,+2,+3])]',
            '[$([CX3,SX3,PX3](=O)[O-,OH])](=O)[O-,OH]',
            '[$([SX4,PX4](=O)(=O)[O-,OH])](=O)(=O)[O-,OH]',
            'c1nn[nH1]n1'
        ],
        "positive charge": [
            'N=[CX3](N)-N',
            '[$([+,+2,+3])&!$(*[-,-2,-3])]',
            '[$([CX3](=N)(-N)[!N])](=N)-N',
            '[$([NX3]([CX4])([CX4,#1])[CX4,#1])&!$([NX3]-*=[!#6])]',
        ]
    }

    smarts_protein = {
        "aromatic ring": [
            "a1aaaa1",
            "a1aaaaa1"
        ],
        "hydrophobicity": [
            '[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,'
            'I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,'
            'I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]',
            '[$([S]~[#6])&!$(S~[!#6])]',
            '[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]',
        ],
        "negative charge": [
            'C(=O)[O-,OH,OX1]',
            '[-,-2,-3,-4]',
        ],
        "positive charge": [
            '[$(C(N)(N)=N)]',
            '[$(n1cc[nH]c1)]',
            '[+,+2,+3,+4]',
        ]
    }

    def __init__(self, file_path):
        self.traj = mdt.load(file_path)
        self.topology = self.traj.topology
        self._coords = puw.quantity(self.traj.xyz, "nanometer")

        self._mol_graph = None
        self._receptor_indices = []

        self._ligand_ids = self.find_ligands()
        self._lig_indices = []
        self._ligand = None
        self._ligand_id = ""

        self._file_path = file_path

    @property
    def ligand_ids(self):
        return self._ligand_ids

    @property
    def ligand(self):
        """ Returns the ligand.

            Returns
            -------
            rdkit.Mol
        """
        return self._ligand

    def set_ligand(self, lig_id):
        """ Set the ligand that will be used to do all computations.

            Parameters
            ----------
            lig_id : str
                Ligand id.
        """
        if lig_id not in self._ligand_ids:
            raise ValueError(
                f"Invalid ligand id: {lig_id}. "
                f"Ligands in this complex:{self._ligand_ids}")
        self._ligand_id = lig_id

    def _update_traj(self, traj):
        """ Update traj, topology and coords attributes.

            Parameters
            ----------
            traj : mdt.Trajectory

        """
        self.traj = traj
        self.topology = traj.topology
        self._coords = puw.quantity(traj.xyz, "nanometer")

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

    def ligand_and_receptor_indices(self):
        """ Get the indices of the ligand with the given id and
            those of the receptor.
        """
        ligand, chain = self._ligand_id.split(":")
        chain_index = self.chain_names.index(chain)
        for atom in self.topology.atoms:
            if atom.residue.name == ligand and atom.residue.chain.index == chain_index:
                self._lig_indices.append(atom.index)
            else:
                self._receptor_indices.append(atom.index)

    def ligand_to_mol(self):
        """ Extract the ligand from the trajectory and create and rdkit mol.

            This will also update the ligand and receptor indices in the complex.
        """
        if len(self._lig_indices) == 0:
            self.ligand_and_receptor_indices()

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

    def remove_ligand(self):
        """ Remove a ligand from the trajectory.
        """
        if len(self._lig_indices) == 0:
            self.ligand_and_receptor_indices()

        self._update_traj(self.traj.atom_slice(self._receptor_indices))
        self._lig_indices.clear()
        self._ligand_ids.remove(self._ligand_id)

    def has_hydrogens(self):
        """ Returns true if the topology contains hydrogens.
        """
        for atom in self.topology.atoms:
            if atom.element.symbol == "H":
                return True
        return False

    def fix_ligand(self, smiles="", add_hydrogens=True):
        """ Correct the ligand bond order and add hydrogens if specified.

            If a smiles is not given, a smiles for the ligand will
            be searched for. In case it is not found this will raise
            an error.

            This method should be called before finding chemical features,
            otherwise they will be incorrect.

            Parameters
            ----------
            smiles : str, optional
                The smiles of the ligand.

            add_hydrogens : bool
                Whether to add hydrogens to the molecule

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

        fixed_lig = Chem.AssignBondOrdersFromTemplate(template, self._ligand)
        if add_hydrogens:
            self._ligand = Chem.AddHs(fixed_lig, addCoords=True)
        else:
            self._ligand = fixed_lig

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

    def add_hydrogens(self):
        """ Add hydrogens to the receptor.

            Necessary to get hydrogen bond protein-ligand interactions.
        """
        modeller = Modeller(
            topology=self.topology.to_openmm(),
            positions=self.traj.openmm_positions(0)
        )
        modeller.addHydrogens()
        self._update_traj(self._modeller_to_trajectory(modeller))

    @staticmethod
    def _mol_to_traj(mol):
        """ Transform a rdkit.Mol to mdtraj.Trajectory

            Parameters
            ----------
            mol : rdkit.Chem.Mol

            Returns
            -------
            traj : mdtraj.Trajectory
        """
        # TODO: convert directly from mol to traj without using file.
        mol_file = tempfile.NamedTemporaryFile(suffix=".pdb")
        Chem.MolToPDBFile(mol, mol_file.name)
        mol_file.seek(0)

        traj = mdt.load(mol_file.name)
        mol_file.close()
        # TODO: converting to traj changes the name of the ligand to UNL
        return traj

    def add_fixed_ligand(self):
        """ Adds the fixed ligand back to the receptor.

            This method should be called after fix_ligand if the original
            topology didn't have hydrogens.
        """
        lig_traj = self._mol_to_traj(self.ligand)

        lig_name = lig_traj.topology.atom(0).residue.name
        n_chains = self.topology.n_chains
        # The ligand will be added to a new chain
        lig_name += ":" + self.chain_names[n_chains]

        coords = np.concatenate((self.traj.xyz, lig_traj.xyz), axis=1)
        topology = self.traj.topology.join(lig_traj.topology)

        self._update_traj(mdt.Trajectory(coords, topology))
        self._ligand_ids = self.find_ligands()
        self.set_ligand(lig_name)
        self.ligand_and_receptor_indices()

    def show(self):
        """ Returns a view of the complex. """
        return show_mdtraj(self.traj)

    def _lig_centroid(self, frame):
        """ Returns the centroid of the ligand at the given frame.

            Parameters
            ----------
            frame : int
                Frame of the trajectory

            Returns
            -------
            puw.Quantity
                Centroid, quantity of shape (3,)
        """
        return np.mean(self._coords[frame, self._lig_indices, :], axis=0)

    def _lig_max_extent(self, centroid, frame):
        """ Returns the maximum extent of the ligand. This is the maximum distance
            from the centroid to any of its atoms.

            Parameters
            ----------
            centroid : puw.Quantity
                Centroid, quantity of shape (3,)

            frame : int
                Frame of the trajectory

            Returns
            -------
            puw.Quantity
                The maximum extent in angstroms, scalar.
       """
        return maths.maximum_distance(
            centroid, self._coords[frame, self._lig_indices, :])

    def binding_site_indices(self, frame):
        """ Obtain the indices of the atoms in the binding site. This corresponds to the
            atoms in the sphere with radius BS_DIST_MAX + ligand maximum extent
            with center in the ligand centroid.

            Parameters
            ----------
            frame : int
                Frame of the trajectory

        """
        self.ligand_and_receptor_indices()
        lig_center = self._lig_centroid(frame)
        lig_extent = self._lig_max_extent(lig_center, frame)
        bs_cutoff = lig_extent + self.BS_DIST_MAX
        distance = np.sqrt(np.sum(np.power(self._coords[0] - lig_center, 2), axis=1))
        return np.where(np.logical_and(distance > lig_extent, distance <= bs_cutoff))[0]

    def ligand_features(self, feat_name, frame):
        """ Returns the centroids and indices of a chemical feature in the ligand.

            Hydrogen bonds cannot be obtained with this method.

            Parameters
            ----------
            feat_name : str
                Name of the feature
            frame : int
                Frame of the trajectory

            Returns
            -------
            centers : list[puw.Quantity]
                List with quantities of shape (3, ).

            indices_list : list[list[int]]
                The indices of the atoms in the topology that correspond
                to the chemical features.

        """
        if self.ligand is None:
            raise exc.NoLigandError

        centers = []
        indices = self.feature_indices(
            self.smarts_ligand[feat_name], self._ligand)
        indices_list = []
        for indices_set in indices:
            feat_coords = self._coords[frame, indices_set, :]
            centers.append(np.mean(feat_coords, axis=0))
            indices_list.append(list(indices_set))

        return centers, indices_list

    def receptor_features(self, feat_name, frame):
        """ Get the centroid and indices of a chemical feature contained
            in the binding site of the receptor.

            Hydrogen bonds cannot be obtained with this method.

            Parameters
            ----------
            feat_name : str
                Name of the feature

            frame : int
                Frame of the trajectory

            Returns
            -------
            centers : list[puw.Quantity]
                List with quantities of shape (3, ).

            indices : list[list[int]]
                The indices of the atoms in the topology that correspond
                to the chemical features.

        """
        if self._mol_graph is None:
            self._mol_graph = Chem.MolFromPDBFile(self._file_path)
            if self._mol_graph is None:
                raise exc.MolGraphError(self._file_path)

        lig_center = self._lig_centroid(frame)
        lig_extent = self._lig_max_extent(lig_center, frame)
        bs_cutoff = lig_extent + self.BS_DIST_MAX

        centers = []
        indices = self.feature_indices(
            self.smarts_ligand[feat_name], self._mol_graph)
        indices_bs = []

        for indices_set in indices:
            feat_coords = self._coords[frame, indices_set, :]
            centroid = np.mean(feat_coords, axis=0)
            # Only keep features in the binding site
            if maths.points_distance(centroid, lig_center) < bs_cutoff:
                centers.append(centroid)
                indices_bs.append(list(indices_set))

        return centers, indices_bs

    @staticmethod
    def feature_indices(feat_def, mol):
        """ Get the indices of the atoms that encompass a chemical
            feature.

            Parameters
            ----------
            feat_def : list[str]
                A smart features definition to find chemical features.

            mol : rdkit.Chem.Mol
                Molecule that will be scanned for the desired features.

            Returns
            -------
            list[tuple[int]]
        """
        feat_indices = []

        for smarts in feat_def:
            pattern = Chem.MolFromSmarts(smarts)
            assert pattern is not None, f"{smarts}"
            all_indices = mol.GetSubstructMatches(pattern)
            for indices in all_indices:
                feat_indices.append(indices)

        return feat_indices

    def hbond_indices(self, frame, criterion="baker"):
        """ Get the indices of atoms involved in hydrogen bonding
            between the ligand and the receptor

            Parameters
            ----------
            frame : int
                Frame of the trajectory

            criterion : str
                The criterion used to compute the hydrogen bonds.

            Returns
            -------
            np.ndarray
                Array of shape (n_hbonds, 3)
        """
        if criterion == "baker":
            # TODO: This function searches for hbonds in all the pl complex
            #   we are only interested in protein-ligand hbonds, we can write
            #   a new function that only searches hbonds between the protein and the ligand.
            return mdt.baker_hubbard(self.traj[frame])
        else:
            raise NotImplementedError

    def hbonds_acceptors(self, hbond_indices, frame):
        """ Returns the coordinates of the hydrogen bond acceptors between the ligand
            and the receptor.

            Parameters
            ----------
            hbond_indices : numpy.ndarray
                Array of shape (n_hbonds, 3)

             frame : int
                Frame of the trajectory

            Returns
            -------
            acceptors : list[puw.Quantity]

        """
        if len(self._lig_indices) == 0:
            self.ligand_and_receptor_indices()

        acceptors = []
        for hbond in hbond_indices:
            don_idx, _, acc_idx = hbond
            if don_idx in self._lig_indices or acc_idx in self._lig_indices:
                acceptors.append(self._coords[frame, acc_idx, :])

        return acceptors

    def hbonds_donors(self, hbond_indices, frame):
        """ Returns the coordinates of the hydrogen bond acceptors between the ligand
           and the receptor.

           Parameters
           ----------
           hbond_indices : numpy.ndarray
               Array of shape (n_hbonds, 3)

            frame : int
               Frame of the trajectory

           Returns
           -------
           acceptors : list[puw.Quantity]

               """
        if len(self._lig_indices) == 0:
            self.ligand_and_receptor_indices()

        donors = []
        for hbond in hbond_indices:
            don_idx, _, acc_idx = hbond
            if don_idx in self._lig_indices or acc_idx in self._lig_indices:
                donors.append(self._coords[frame, don_idx, :])

        return donors

    def interactions_view(self, bsite_indices, feats=None):
        """ Returns a view of the binding site and the ligand,
            optionally with the specified chemical features highlighted
            as spheres.

            Parameters
            ----------
            bsite_indices : list[int] or numpy.ndarray
                A list of indices that define the binding site.

            feats : list[dict[str, list[puw.Quantity]]
                A list of dictionaries containing a list of the
                coordinates of each feature type.

            Returns
            -------
            view : nglview.NGLWidget
        """
        palette = {
            'positive charge': '#3498DB',  # Blue
            'negative charge': '#884EA0',  # Purple
            'hb acceptor': '#B03A2E',  # Red
            'hb donor': '#17A589',  # Green
            'hydrophobicity': '#F5B041',  # Orange
            'aromatic ring': '#F1C40F',  # Yellow
        }

        bsite_traj = self.traj.atom_slice(bsite_indices)
        view = show_mdtraj(bsite_traj)

        for feats_dict in feats:
            for feat_name, coord_list in feats_dict.items():
                for coord in coord_list:
                    value = puw.get_value(coord, "angstroms").tolist()
                    view.shape.add_sphere(
                        value, to_rgb(palette[feat_name]), 1.0, feat_name
                    )
                    n_components = len(view._ngl_component_ids)
                    view.update_representation(
                        component=n_components, repr_index=0, opacity=0.5)

        return view

    def prepare(self, lig_id, smiles="", add_hydrogens=True):
        """ Fix the ligand bond order, and optionally add hydrogens to
            both the ligand and the receptor.

            Adding hydrogens is necessary to find hydrogen bonds.

            If a smiles is not given, a smiles for the ligand will
            be searched for. In case it is not found this will raise
            an error.

            Parameters
            ----------
            lig_id: str
                Id of the ligand that will be fixed.

            smiles : str, optional
                The smiles of the ligand.

            add_hydrogens: bool, optional.
                Whether to add hydrogens to the ligand and the receptor.

        """
        pass
