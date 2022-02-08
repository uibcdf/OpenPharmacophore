# OpenPharmacophore
from openpharmacophore._private_tools.exceptions import InvalidFileFormat, NoLigandsError, OpenPharmacophoreTypeError
from openpharmacophore.pharmacophore.pharmacophoric_point import UniquePharmacophoricPoint
from openpharmacophore import StructuredBasedPharmacophore
from openpharmacophore import Pharmacophore
from openpharmacophore.utils.conformers import conformer_energy
from openpharmacophore.pharmacophore.color_palettes import get_color_from_palette_for_feature
# Third Party
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.lib.util import NamedStream
import mdtraj as mdt
import numpy as np
import pandas as pd
import pyunitwizard as puw
from rdkit.Chem.Draw import rdMolDraw2D
from tqdm.auto import tqdm
# Standard Library
from collections import defaultdict
import copy
import bisect
from io import StringIO
import tempfile
from typing import List, Tuple, Optional

class Dynophore():
    """ Class to store and compute dynamic pharmacophores

    Parameters
    ----------

    trajectory : str or mdtraj.trajectory or MDAnalysis.universe
        A str with the file path containing the trajectory, an mdtraj trajectory object, 
        or an MDAnalysis universe.

    Attributes
    ----------

    pharmacophores : list of openpharmacophore.StructuredBasedPharmacophore
        List with pharmacophores for each relevant frame in the trajectory. 

    pharmacophore_indices : list of int
        Indices of the frame of the trajectory from which the pharmacophores were extracted.
        The index of each element of the list corresponds to the one in pharmacophores list.

    n_pharmacophores : int
        Number of different pharmacophores in the trajectory.

    """
    def __init__(self, trajectory):
        self.pharmacophores = []
        self.pharmacophore_indices = []
        self.n_pharmacophores = 0
        self.unique_pharmacophoric_points = []

        # TODO: Load other types of file, including using a topology and tajectory
        if isinstance(trajectory, str):
            self._trajectory = self._load_trajectory_file(trajectory)
        elif isinstance(trajectory, mdt.Trajectory):
            self._trajectory_type = "mdt"
            self._trajectory = trajectory
            self._n_frames = self._trajectory.n_frames
        elif isinstance(trajectory, mda.Universe):
            self._trajectory_type = "mda"
            self._trajectory = trajectory
            self._n_frames = trajectory.trajectory.n_frames
        else:
            raise TypeError("Trajectory must be of type string, mdtraj.Trajectory or MdAnalysis.Universe")
        
        
        self._saved_ligand = False
        self._averaged_coords = False

    def common_hits_approach(self, frame_list=None):
        """ Get a list of pharmacophore models from a trajectory using the common hits approach
            method.

            Notes
            -----

            This method is based on obtaining a list of representative pharmacophore models from a 
            trajectory and then validate and score them using virtual screening. The best performant
            pharmacophore models are then returned.

            References
            ----------

            [1] Wieder, Marcus, Arthur Garon, Ugo Perricone, Stefan Boresch, Thomas Seidel, Anna Maria Almerico, 
            and Thierry Langer. "Common hits approach: combining pharmacophore modeling and molecular dynamics 
            simulations." Journal of chemical information and modeling 57, no. 2 (2017): 365-385      

        """
        if frame_list is None:
            frame_list = list(range(0, self._n_frames))

        self.pharmacophores_from_frames(frame_list, load_ligand=True)
        self._get_unique_pharmacophoric_points(avg_coordinates=False)
        rpms = self.representative_pharmacophore_models()

        pass

    def draw(self, file_name: str, img_size: Tuple[int, int] = (500,500), 
            legend: str = "", freq_threshold: float = 0.2) -> None:
        """ Draw a 2d representation of the dynamic pharmacophore. This is a drawing of the
            ligand with the pharmacophoric features highlighted and the frequency if each
            one. 

            Parameters
            ----------
            file_name : str
                Name or path og the file where the drawing will be saved. Must be a png file.

            img_size : 2-tuple of int, optional 
                The size of the image (default=(500,500))

            legend : str, optional
                Image legend.

            freq_threshold : double , optional
                The minimun frequency of a pharmacophoric point to be drawn. Number
                between 0.0 and 1.0 (default=0.2). 
        """
        if freq_threshold < 0.0 or freq_threshold > 1.0:
            raise ValueError("Freqency threshold must be a value between 0 and 1")    

        if not file_name.endswith(".png"):
            raise InvalidFileFormat("File must be a png.")

        # Extract a ligand
        if self.pharmacophores[0].ligand is None:
            raise NoLigandsError("Ligand could not be extracted")
        ligand = copy.deepcopy(self.pharmacophores[0].ligand)
        ligand.RemoveAllConformers()

        atoms = []
        bond_colors = {}
        atom_highlights = defaultdict(list)
        highlight_radius = {}

        for up in self.unique_pharmacophoric_points:
            
            if up.frequency < freq_threshold:
                continue

            indices = up.atom_indices
            update_freq = True
            for idx in indices:

                # If an atom has more than one feature keep higher frequency value
                if idx in atoms:
                    if ligand.GetAtomWithIdx(idx).HasProp("atomNote"):
                        freq = int(ligand.GetAtomWithIdx(idx).GetProp("atomNote")[2:])
                        if freq > up.frequency:
                            update_freq = False

                atoms.append(idx)
                if "hydrophobicity" in up.feature_name:
                    feat_name = "hydrophobicity"
                else:
                    feat_name = " ".join(up.feature_name.split()[0:2])
                
                atom_highlights[idx].append(get_color_from_palette_for_feature(feat_name))
                highlight_radius[idx] = 0.6

                # Draw aromatic rings bonds
                if up.short_name == "R":
                    for neighbor in ligand.GetAtomWithIdx(idx).GetNeighbors():
                        nbr_idx = neighbor.GetIdx()
                        if nbr_idx not in indices:
                            continue
                        bond = ligand.GetBondBetweenAtoms(idx, nbr_idx).GetIdx()
                        bond_colors[bond] = [get_color_from_palette_for_feature("aromatic ring")]
            
            if update_freq:
                frequency = int(up.frequency * 100)
                ligand.GetAtomWithIdx(idx).SetProp("atomNote", f"f={frequency}")

        drawing = rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])
        drawing.DrawMoleculeWithHighlights(ligand, legend, dict(atom_highlights), bond_colors, highlight_radius, {})
        drawing.FinishDrawing()
        drawing.WriteDrawingText(file_name)

    def first_and_last_pharmacophore(self) -> None:
        """ Derive a pharmacophore model for the first and last frames of a trajectory.

            References
            ----------
            [1] Wieder, Marcus, Ugo Perricone, Thomas Seidel, Stefan Boresch, and Thierry Langer. 
            "Comparing pharmacophore models derived from crystal structures and from molecular 
            dynamics simulations." Monatshefte für Chemie-Chemical Monthly 147, no. 3 (2016): 
            553-563.
        """
        if self._trajectory_type == "mdt":
            get_pharmacophore = self._pharmacophore_from_mdtraj
        elif self._trajectory_type == "mda":
            get_pharmacophore = self._pharmacohore_from_mdanalysis

        initial_pharmacophore = get_pharmacophore(0, True, True)
        end_pharmacophore = get_pharmacophore(-1, True, True)
        last_frame_index = self._trajectory.n_frames
        self.pharmacophores = [
            initial_pharmacophore,
            end_pharmacophore
        ]
        self.pharmacophore_indices = [0, last_frame_index]
        self.n_pharmacophores = 2

    def pharmacophore_by_frequency(self, threshold: float) -> Pharmacophore:
        """ Derive a unique pharmacophore model with the pharmacophoric points
            that have a frequency >= to threshold.

            Parameters
            ---------
            threshold : float
                The value of frequency from which points are considered part of
                the pharmacophore model. Must be a value between 0 and 1-

            Returns
            -------
            openpharmcophore.Pharmacophore
                Pharmacophore model with the unique pharmacophoric points.

            References
            ----------
            [1] Wieder, Marcus, Ugo Perricone, Thomas Seidel, and Thierry Langer. "Pharmacophore models 
            derived from molecular dynamics simulations of protein-ligand complexes: A case study." 
            Natural product communications 11, no. 10 (2016): 1934578X1601101019.
        """
        if threshold < 0 or threshold > 1:
            raise ValueError("Threshold must be a number between 0 and 1")
        
        if len(self.unique_pharmacophoric_points) == 0:
            self._get_unique_pharmacophoric_points(avg_coordinates=True)
            
        points = [p for p in self.unique_pharmacophoric_points if p.frequency >= threshold]
        return Pharmacophore(points)

    def pharmacophore_from_unique_points(self, unique_points: List[str]) -> Pharmacophore:
        """ Get a pharmacophore which consists of the passed unique pharmacophoric
            points.

            Parameters
            ----------
            unique_points: list of str
                List with the name of the unique pharmacophoric points.

            Returns
            -------
            openpharmcophore.Pharmacophore
                Pharmacophore model with the specified points.
        """
        if len(self.unique_pharmacophoric_points) == 0 or not self._averaged_coords:
            self._get_unique_pharmacophoric_points(avg_coordinates=True)
        points = [point for point in self.unique_pharmacophoric_points if point.feature_name in unique_points]
        return Pharmacophore(pharmacophoric_points=points)

    def pharmacophores_from_frames(self, frames: List[int], load_ligand: bool = True) -> None:
        """ Get pharmacophores for the specified frames in a trajectory

            Parameters
            ----------
            frames : list of int
                Indices of the frames for which pharmacophores will be derived.

        """
        if self._trajectory_type == "mdt":
            get_pharmacophore = self._pharmacophore_from_mdtraj
        elif self._trajectory_type == "mda":
            get_pharmacophore = self._pharmacohore_from_mdanalysis
        
        self.pharmacophores.clear()
        self.pharmacophore_indices.clear()
        for ii in tqdm(frames):
            self.pharmacophores.append(get_pharmacophore(ii, load_ligand=load_ligand))
            self.pharmacophore_indices.append(ii)
        self.n_pharmacophores = len(self.pharmacophores)
    
    def pharmacophoric_point_frequency(self) -> pd.DataFrame:
        """ Get a dataframe with all unique pharmacophoric points and its frequency.

            Returns
            -------
            pandas.DataFrame
                Dataframe with the following columns: feature name, frequency and atom
                indices.
        """
        if len(self.unique_pharmacophoric_points) == 0 or not self._averaged_coords:
            self._get_unique_pharmacophoric_points(avg_coordinates=True)
        
        names = []
        frequencies = []
        indices = []
        for point in self.unique_pharmacophoric_points:
            names.append(point.feature_name)
            frequencies.append(point.frequency)
            indices.append(point.atom_indices)

        frequency = pd.DataFrame().from_dict({
            "Feature Name": names,
            "Frequency": frequencies,
            "Atoms Indices": indices
        })
        frequency.sort_values(by=["Frequency"], ascending=False, inplace=True)
        frequency.reset_index(inplace=True)
        frequency.drop(columns=["index"], inplace=True)
        return frequency

    def point_frequency_plot(self, threshold: float = 0.0, n_bins: int = 10, 
            ax: Optional[plt.Axes] = None):
        """ Plot of pharmacophoric points frequency vs time. 
        
            Each pharmacophoric point will appear as a different line in the plot.

            Parameters
            ----------
            threshold : double, default=0.0
                The value of overall frequency from which points will form part of the 
                plot. If there are a lot of points with really low frequency, setting
                the threshold value can help with visualization.

            n_bins : int, default=10
                Number of bins to discretize the timesteps.            

            ax : matplotlib.axes._subplots.AxesSubplot, optional.
                An axes object where the plot will be drawn.
        """
        if len(self.unique_pharmacophoric_points) == 0 or not self._averaged_coords:
            self._get_unique_pharmacophoric_points(avg_coordinates=True)

        if threshold < 0 or threshold > 1:
            raise ValueError("Threshold must be a number between 0 and 1")

        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 7))
        n_timesteps = self._n_frames
        bins = np.arange(0, n_timesteps + 1, n_timesteps/n_bins)

        for point in self.unique_pharmacophoric_points:
            if point.frequency < threshold:
                continue
            point_timesteps = np.array(point.timesteps)
            discretized_timesteps = np.digitize(point_timesteps, bins)

            counts = np.zeros_like(bins)

            for i in range(bins.shape[0]):
                c = np.count_nonzero(discretized_timesteps == i)
                counts[i] = c
            
            ax.plot(bins, counts, label=point.feature_name)

        ax.legend()
        ax.set_xlabel("Timesteps")
        ax.set_ylabel("Count")
        plt.show()

        return ax
    
    def representative_pharmacophore_models(self) -> List[StructuredBasedPharmacophore]:
        """ Get all representative pharmacophore models (RPM) in a trajectory. 
            
            RPMs are pharmacophore models that have the same pharmacophoric points, 

            Returns
            -------
            rpms : list of openpharmacophore.StructuredBasedPharmacophore
                The representative pharmacophore models

            Note
            -----
            Pharmacophoric points are considered equal based only  on feature type and the atoms to 
            which this points belong to. Coordinates are not taken into account.

            The coordinates of the pharmacophoric points are those that belong to the median energy of
            the ligand.

            References
            ----------
            [1] Wieder, Marcus, Arthur Garon, Ugo Perricone, Stefan Boresch, Thomas Seidel, Anna Maria Almerico, 
            and Thierry Langer. "Common hits approach: combining pharmacophore modeling and molecular dynamics 
            simulations." Journal of chemical information and modeling 57, no. 2 (2017): 365-385   
            
        """
        if len(self.unique_pharmacophoric_points) == 0 or self._averaged_coords:
            self._get_unique_pharmacophoric_points(avg_coordinates=False)
            self._averaged_coords = False
        
        rpms_indices = self._get_rpms_indices()
        
        return self._pharmacophores_from_ligand_median_energy(rpms_indices)

    def _get_rpms_indices(self) -> List[List[int]]:
        """ Get the indices of the representative pharmacophore models.
        
            If an empty list is returned it means that all pharmacophore models in the trajectory are different.
        
            Returns
            --------
            rpms_indices : list of list of int
                A list where each sublist contains the indices of each representative pharmacophore
                model. This indices correspond to the attribute pharmacophores of the Dynophore
                class.
        """
        # Compute a matrix where each row represents a feature vector of a pharmacophore
        n_pharmacophores = self.n_pharmacophores
        n_features = len(self.unique_pharmacophoric_points) 
        feature_matrix = np.zeros((n_pharmacophores, n_features), dtype=np.int32)
        for ii, pharmacophore in enumerate(self.pharmacophores):
            for point in pharmacophore:
                for jj, unique_point in enumerate(self.unique_pharmacophoric_points):
                    if point.is_equal(unique_point):
                        feature_matrix[ii, jj] = 1
                        break
        
        # Find similar pharmacophores in the matrix
        rpms_indices = []
        skip = []
        for ii in range(n_pharmacophores):
            rpm = [ii]
            for jj in range(ii + 1, n_pharmacophores):
                if jj in skip:
                    continue
                if np.all(feature_matrix[ii, :] == feature_matrix[jj, :]):
                    rpm.append(jj)
                    skip.append(jj)
            # Keep only models that have a frequency higher than 2
            if len(rpm) > 2:
                rpms_indices.append(rpm)
        
        
        return rpms_indices
       
    def _pharmacophores_from_ligand_median_energy(self, rpms_indices)-> List[List[int]]:
        """ Get the representative pharmacophore models that correspond to the pharmacophore
            with ligand median energy.

            Parameters
            ----------
            rpms_indices : list of list of int
                A list where each sublist contains the indices of each representative pharmacophore
                model. This indices correspond to the attribute pharmacophores of the Dynophore
                class.
            
            Returns
            -------
            rpms : list of openpharmacophore.StructuredBasedPharmacophore
                The representative pharmacophore models
        """
        rpms = []
        for indices in rpms_indices:
            energies = []
            for index in indices:
                energy = (conformer_energy(self.pharmacophores[index].ligand), index)
                bisect.insort(energies, energy)
            # Take the pharmacophore with median energy
            median_energy_index = energies[int(len(energies) / 2)][1]
            rpms.append(self.pharmacophores[median_energy_index])
        
        return rpms

            
    def _load_trajectory_file(self, file_name: str) -> mdt.Trajectory:
        """ Load a trajectory file from a MD simulation

            Parameters
            ----------
            file_name : str
                Name of the file containing the trajectory.

            Returns
            -------
            traj : 
                The trajectory object.  
        """
        if file_name.endswith("h5"):
            traj = mdt.load(file_name)
            self._trajectory_type = "mdt"
        else:
            raise NotImplementedError

        return traj
    
    def _get_unique_pharmacophoric_points(self, avg_coordinates: bool = True) -> None:
        """ Get all unique pharmacophoric points across all the pharmacophore models 
            derived from the trajectory. 

            Parameters
            ----------
            avg_coordinates : bool
                Whether to average the coordinates of the pharmacophoric points.
            
            Notes
            -----
            Two points are considered equal if they have the same feature type and
            are associated with the same atom in the ligand.
        """
        if avg_coordinates:
            self._averaged_coords = True

        if self.n_pharmacophores == 0:
            self.pharmacophores_from_frames(list(range(0, self._n_frames)))
        
        all_points = []
        for ii, pharmacophore in enumerate(self.pharmacophores):
            for pharmacophoric_point in pharmacophore:
                pharmacophoric_point.pharmacophore_index = ii
                all_points.append(pharmacophoric_point)
        
        self.unique_pharmacophoric_points.clear()
        # Get all unique parmacophoric points while also updating the count, 
        # timesteps where they appear and calculating the average centroid.
        for point in all_points:
            is_unique = True
            for unique_p in self.unique_pharmacophoric_points:
                if point.is_equal(unique_p):
                    timestep = point.pharmacophore_index
                    if not timestep in unique_p.timesteps:
                        unique_p.timesteps.append(timestep)
                        unique_p.count += 1
                    if avg_coordinates:
                        unique_p.center += point.center
                    is_unique = False
                    break
            if is_unique:
                self.unique_pharmacophoric_points.append(UniquePharmacophoricPoint(point, point.pharmacophore_index))
        
        names = []
        for point in self.unique_pharmacophoric_points:
            if avg_coordinates:
                # Normalize centroid
                point.center /= point.count 
            point.frequency = point.count / self.n_pharmacophores
            # Get a unique name for each point
            feat_num = 1
            full_name = point.feature_name + " " + str(feat_num)
            if full_name not in names:
                names.append(full_name)
                point.feature_name = full_name
            else:
                while True:
                    feat_num += 1
                    full_name = point.feature_name + " " + str(feat_num)
                    if full_name not in names:
                        names.append(full_name)
                        point.feature_name = full_name
                        break

    def _pharmacophore_from_mdtraj(self, frame_num: int, load_mol_system: bool=False, 
                                load_ligand: bool=False) -> StructuredBasedPharmacophore:
        """ Derive a pharmacophore for a single frame of an mdtraj Trajectory object.

            Parameters
            ----------
            frame_num : int
                The index number of the frame from which the pharmacophore will be derived.
            
            load_mol_system : bool, default=False
                If true the receptor will be stored in the pharmacophore object.
            
            load_ligand : bool, default=False
                If true the ligand will be stored in the pharmacophore object.
        """
        # mdtraj trajectories cannot be passed to SringIO objects nor saved as string. So with this
        # method, temporary pdb files will be created that can be read by the StructuredBasedPharmacophore 
        # class.
        if not isinstance(frame_num, int):
            raise OpenPharmacophoreTypeError("Frame number must be an integer")
        frame = self._trajectory[frame_num]
        
        with tempfile.NamedTemporaryFile() as original_file:
            frame.save_pdb(original_file.name)
            original_file.seek(0) 
            lines_original = original_file.readlines()
                
        # The pdb mdtraj generates needs to be edited so that pybel can read it.
        # The third line that contains "MODEL" needs to be removed for the structured 
        # based pharmacophore to work.
        with tempfile.NamedTemporaryFile() as modified_file:
            for line in lines_original:
                if not line.startswith(b'MODEL'):
                    modified_file.write(line)
            modified_file.truncate()
            modified_file.seek(0)
            pharmacophore = StructuredBasedPharmacophore.from_pdb(modified_file, 
                    radius=1.0, ligand_id=None, hydrophobics="plip", 
                load_mol_system=load_mol_system, load_ligand=load_ligand)
        
        return pharmacophore
    
    def _pharmacohore_from_mdanalysis(self, frame_num: int, load_mol_system: bool = False, 
                                    load_ligand: bool = False) -> StructuredBasedPharmacophore:
        """ Derive a pharmacophore for a single frame of an MdAnalysis Universe object.

            Parameters
            ----------
            frame_num : int
                The index number of the frame from which the pharmacophore will be derived.
            
            load_mol_system: bool, default=False
                If true the receptor will be stored in the pharmacophore object.
            
            load_ligand: bool, default=False
                If true the ligand will be stored in the pharmacophore object.
        """
        if not isinstance(frame_num, int):
            raise OpenPharmacophoreTypeError("Frame number must be an integer")
        stream = StringIO()
        pdb_stream = NamedStream(stream, "output.pdb")
        atoms = self._trajectory.select_atoms("all")
        atoms.write(pdb_stream, frames=self._trajectory.trajectory[[frame_num]])
        pharmacophore = StructuredBasedPharmacophore.from_pdb(pdb_stream, 
            radius=1.0, ligand_id=None, hydrophobics="plip", 
            load_mol_system=load_mol_system, load_ligand=load_ligand)
        
        return pharmacophore

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(n_pharmacophores={self.n_pharmacophores}; n_frames={self._n_frames})"

    
    
