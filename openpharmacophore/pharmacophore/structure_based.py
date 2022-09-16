from .pharmacophore import Pharmacophore
from .pharmacophoric_point import PharmacophoricPoint
from .rdkit_pharmacophore import rdkit_pharmacophore
from ..io import (json_pharmacophoric_elements, ligandscout_xml_tree,
                  mol2_file_info, ph4_string)
from ..io import (load_json_pharmacophore, load_mol2_pharmacophoric_points,
                  pharmacophoric_points_from_ph4_file, read_ligandscout)
from .._private_tools.exceptions import InvalidFileFormat, PDBFetchError
from plip.structure.preparation import PDBComplex
import numpy as np
import nglview as nv
import pyunitwizard as puw
import json
import re
import requests


class StructureBasedPharmacophore(Pharmacophore):
    """ Class to store, and extract pharmacophores from protein-ligand complexes.

        The pharmacophores can be extracted from a pdb file or from a molecular
        dynamics simulation.

    """

    def __init__(self, receptor):
        # Pharmacophores will be stored as a list of pharmacophoric points.
        # A list for each pharmacophore
        self._pharmacophores = []
        self._pharmacophores_frames = []  # Contains the frame to which each pharmacophore belongs
        self._topology = None
        self._coordinates = None
        self._num_frames = 0

        if isinstance(receptor, str):
            if self._is_pdb_id(receptor):
                pdb_str = self._fetch_pdb(receptor)
                self.extract(pdb_str, True)
                self._num_frames += 1
            elif receptor.endswith(".pdb"):
                self.extract(receptor, False)
                self._num_frames += 1
            else:
                # Load from a trajectory
                raise NotImplementedError

    @property
    def num_frames(self):
        return self._num_frames

    @staticmethod
    def _is_pdb_id(receptor):
        """ Check if the receptor is a PDB id.

            Parameters
            ----------
            receptor: str
                The receptor should be a string

            Returns
            -------
            bool
                Whether the receptor is a pdb id.
        """
        if len(receptor) == 4:
            pattern = re.compile('[0-9][a-zA-Z_0-9]{3}')
            if pattern.match(receptor):
                return True
        return False

    @staticmethod
    def _fetch_pdb(pdb_id):
        """ Fetch a PDB with the given id.
        """
        url = f'http://files.rcsb.org/download/{pdb_id}.pdb'
        res = requests.get(url, allow_redirects=True)

        if res.status_code != 200:
            raise PDBFetchError(pdb_id, url)

        return res.content.decode()

    def from_file(self, file_name):
        """ Load a pharmacophore from a file.

          Parameters
          ---------
          file_name : str
              Name of the file containing the pharmacophore

        """
        if file_name.endswith(".json"):
            self._pharmacophores.append(load_json_pharmacophore(file_name)[0])
            self._num_frames = 1
            self._pharmacophores_frames.append(0)
        elif file_name.endswith(".mol2"):
            # Loads all pharmacophores contained in the mol2 file. It assumes that each pharmacophore
            # corresponds to a different timestep
            self._pharmacophores = load_mol2_pharmacophoric_points(file_name)
            self._pharmacophores_frames = list(range(0, len(self._pharmacophores)))
            self._num_frames = len(self._pharmacophores)
        elif file_name.endswith(".pml"):
            self._pharmacophores.append(read_ligandscout(file_name))
            self._num_frames = 1
            self._pharmacophores_frames.append(0)
        elif file_name.endswith(".ph4"):
            self._pharmacophores.append(pharmacophoric_points_from_ph4_file(file_name))
            self._num_frames = 1
            self._pharmacophores_frames.append(0)
        else:
            raise InvalidFileFormat(file_name.split(".")[-1])

    def add_frame(self):
        """ Add a new frame to the pharmacophore. """
        self._pharmacophores.append([])
        self._pharmacophores_frames.append(self._num_frames)
        self._num_frames += 1

    def add_points_to_frame(self, point_list, frame):
        """ Add pharmacophoric points from a list to a frame. """
        for point in point_list:
            self._pharmacophores[frame].append(point)

    def add_point(self, point, frame):
        """ Add a pharmacophoric point to a pharmacophore in a specific frame."""
        self._pharmacophores[frame].append(point)

    def remove_point(self, index, frame):
        """ Removes a pharmacophoric point from the pharmacophore at the given frame."""
        self._pharmacophores[frame].pop(index)

    def remove_picked_point(self, view):
        pass

    def edit_picked_point(self, view):
        pass

    def add_point_in_picked_location(self, view):
        pass

    def add_to_view(self, view, frame=None):
        """ Add pharmacophore(s) to a ngl view.
        """
        if frame is None or isinstance(frame, list):
            raise NotImplementedError
        else:
            for point in self[frame]:
                point.add_to_ngl_view(view)

    def show(self, frame=None):
        """ Shows a 3D representation of the pharmacophore model. """
        view = nv.NGLWidget()
        self.add_to_view(view, frame)

    def to_json(self, file_name, frame):
        """ Save pharmacophore(s) to a json file.
        """
        data = json_pharmacophoric_elements(self[frame])
        with open(file_name, "w") as fp:
            json.dump(data, fp)

    def to_ligand_scout(self, file_name, frame):
        """ Save a pharmacophore at a given frame to ligand scout format (pml).
        """
        xml_tree = ligandscout_xml_tree(self[frame])
        xml_tree.write(file_name, encoding="UTF-8", xml_declaration=True)

    def to_moe(self, file_name, frame):
        """ Save a pharmacophore at a given frame to moe format (ph4).
        """
        pharmacophore_str = ph4_string(self[frame])
        with open(file_name, "w") as fp:
            fp.write(pharmacophore_str)

    def to_mol2(self, file_name, frame=None):
        """ Save pharmacophore(s) to mol2 file.
        """
        # TODO: save multiple pharmacophores
        if frame is None or isinstance(frame, list):
            raise NotImplementedError
        pharmacophore_data = mol2_file_info([self[frame]])
        with open(file_name, "w") as fp:
            fp.writelines(pharmacophore_data[0])

    def to_rdkit(self, frame):
        """ Transform a pharmacophore at a given frame to a rdkit pharmacophore.
        """
        return rdkit_pharmacophore(self[frame])

    def extract(self, receptor, frames=None):
        """ Extract pharmacophore(s) from the receptor.
        """
        pass

    @staticmethod
    def _protein_ligand_interactions(receptor, as_string=False):
        """ Returns a dictionary with the protein-ligand interactions.

            Parameters
            ----------
            receptor : str
                File or string containing the protein-ligand complex.

            as_string : bool
                Whether the receptor is a string or a file path.

            Returns
            -------
            ligands : dict[str, plip.PLInteraction]
                The interactions of each ligand in the receptor
        """
        protein = PDBComplex()
        protein.load_pdb(receptor, as_string=as_string)
        protein.analyze()
        # Remove ions and very small molecules from the interactions dictionary
        interactions = protein.interaction_sets
        for ligand in protein.ligands:
            if ligand.type != "SMALLMOLECULE" or len(ligand.atomorder) <= 5:
                ligand_id = ":".join([ligand.hetid, ligand.chain, str(ligand.position)])
                interactions.pop(ligand_id)

        return interactions

    @staticmethod
    def _points_from_interactions(interactions, ligand, radius):
        """ Transform interactions to a list of pharmacophoric points.

            Parameters
            ----------
            interactions : dict[str, plip.PLInteraction]
                Interaction data for each ligand in the receptor.

            Returns
            -------
            points : list[PharmacophoricPoint]
                The pharmacophoric points

        """
        interactions_list = interactions[ligand].all_itypes
        points = []
        radius = puw.quantity(radius, "angstroms")

        for interaction in interactions_list:
            interaction_name = type(interaction).__name__

            if interaction_name == "pistack":
                ligand_center = np.array(interaction.ligandring.center)
                protein_center = np.array(interaction.proteinring.center)
                direction = protein_center - ligand_center
                atom_indices = {atom.idx for atom in interaction.ligandring.atoms}
                aromatic = PharmacophoricPoint(
                    feat_type="aromatic ring",
                    center=puw.quantity(ligand_center, "angstroms"),
                    radius=radius,
                    direction=direction,
                    atom_indices=atom_indices
                )
                points.append(aromatic)

            elif interaction_name == "hydroph_interaction":
                center = puw.quantity(interaction.ligatom.coords, "angstroms")
                atom_inx = {interaction.ligatom.idx}
                hydrophobic = PharmacophoricPoint(
                    feat_type="hydrophobicity",
                    center=center,
                    radius=radius,
                    direction=None,
                    atom_indices=atom_inx
                )
                points.append(hydrophobic)

            elif interaction_name == "saltbridge":
                if interaction.protispos:
                    # The ligand has a negative charge
                    center = puw.quantity(interaction.negative.center, "angstroms")
                    atom_indices = {atom.idx for atom in interaction.negative.atoms}
                    charge_sphere = PharmacophoricPoint(
                        feat_type="negative charge",
                        center=center,
                        radius=radius,
                        atom_indices=atom_indices
                    )
                else:
                    # The ligand has a positive charge
                    center = puw.quantity(interaction.positive.center, "angstroms")
                    atom_indices = {atom.idx for atom in interaction.positive.atoms}
                    charge_sphere = PharmacophoricPoint(
                        feat_type="positive charge",
                        center=center,
                        radius=radius,
                        atom_indices=atom_indices
                    )
                points.append(charge_sphere)

            elif interaction_name == "hbond":
                if interaction.protisdon:
                    # The ligand has an acceptor atom
                    ligand_acceptor_center = np.array(interaction.a.coords)
                    ligand_acceptor_inx = {interaction.a.idx}
                    protein_donor_center = np.array(interaction.d.coords)
                    direction = ligand_acceptor_center - protein_donor_center
                    acceptor = PharmacophoricPoint(
                        feat_type="hb acceptor",
                        center=puw.quantity(ligand_acceptor_center, "angstroms"),
                        radius=radius,
                        direction=direction,
                        atom_indices=ligand_acceptor_inx
                    )
                    points.append(acceptor)
                else:
                    # The ligand has a donor atom
                    ligand_donor_center = np.array(interaction.d.coords)
                    ligand_donor_inx = {interaction.d.idx}
                    protein_acceptor_center = np.array(interaction.a.coords)
                    direction = protein_acceptor_center - ligand_donor_center
                    donor = PharmacophoricPoint(
                        feat_type="hb donor",
                        center=puw.quantity(ligand_donor_center, "angstroms"),
                        radius=radius,
                        direction=direction,
                        atom_indices=ligand_donor_inx
                    )
                    points.append(donor)
            else:
                continue

        return points

    def __len__(self):
        return len(self._pharmacophores)

    def __getitem__(self, frame):
        return self._pharmacophores[frame]
