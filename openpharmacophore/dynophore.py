from openpharmacophore.structured_based import StructuredBasedPharmacophore
from openpharmacophore.pharmacophore import Pharmacophore
from openpharmacophore.utils.random_string import random_string
import mdtraj as md
import os

class Dynophore():
    """ Class to store and compute dynamic pharmacophores

    Parameters
    ----------

    trajectory : 
        A file containing the trajectory or a trajectory object.

    Attributes
    ----------

    pharmacophores : dict
        Dictionary with pharmacophores for each relevant frame in the
        trajectory. 

    n_pharmacophores : int
        Number of different pharmacophores in the trajectory.

    """
    def __init__(self, trajectory, method="frequency"):
        self.pharmacophores = {}
        self.n_pharmacophores = 0
        self.method = method.lower()

        if isinstance(trajectory, str):
            self.trajectory = self._load_trajectory_file(trajectory)
        else:
            self.trajectory = trajectory

    def get_pharmacophores(self):
        if self.method == "last":
            initial_pharmacophore = self._pharmacophore_from_frame(0)
            end_pharmacophore = self._pharmacophore_from_frame(-1)
            last_frame_index = self.trajectory.n_frames
            self.pharmacophores = {
                "frame-0": initial_pharmacophore,
                f"frame-{last_frame_index}": end_pharmacophore
            }
        elif self.method == "frequency":
            self._pharmacophores_by_frequency()
        else:
            raise NotImplementedError 

    def _load_trajectory_file(self, file_name):
        """ Load a trajectory file from a MD simulation

            Parameters
            ----------
            file_name: str
                Name of the file containing the trajectory.

            Returns
            -------
            traj: 
                The trajectory object.  
        """
        if file_name.endswith("h5"):
            traj = md.load(file_name)
        else:
            raise NotImplementedError

        return traj
    
    def _pharmacophore_from_frame(self, frame_num):
        """ Derive a pharmacophore for a single frame of the trajectory
        """
        if not isinstance(frame_num, int):
            raise ValueError("Frame number must be an integer")
        temp_filename = "./temp" + random_string(10) + ".pdb"
        frame = self.trajectory[frame_num]
        frame.save_pdb(temp_filename)

        # The pdb mdtraj generates needs to be edited so that pybel can read it.
        # The third line that contains "MODEL" needs to be removed for the structured 
        # based pharmacophore to work.
        with open(temp_filename, "r+") as f:
            d = f.readlines()
            f.seek(0)
            for i in d:
                if  not i.startswith("MODEL"):
                    f.write(i)
            f.truncate()
        pharmacophore = StructuredBasedPharmacophore.from_pdb(temp_filename, 
            radius=1.0, ligand_id=None, hydrophobics="plip", load_mol_system=True)
        
        os.remove(temp_filename)
        return pharmacophore
    
    def _pharmacophores_by_frequency(self):
        pharmacophores = []
        for i in range(self.trajectory.n_frames):
            pharmacophores.append(self._pharmacophore_from_frame(i))

    # Methods for obtaining dynophores:
    # Last-frame: consider only the pharmacophores form the first and last frame of the MD trajectory.
    #             See:  Wieder, Marcus, Ugo Perricone, Thomas Seidel, Stefan Boresch, and Thierry Langer. 
    #                   "Comparing pharmacophore models derived from crystal structures and from molecular 
    #                   dynamics simulations." Monatshefte für Chemie-Chemical Monthly 147, no. 3 (2016): 
    #                   553-563.
    #
    # Frequency: 
    #       1.  Derive a pharmacophore model for each frame.
    #       2.  Count the frequency of each feature. For the frequency we don't take the coordinates
    #           into account. Just the feature type and the atoms that are part of that feature.
    #       3.  Calculate frequency in terms of percentage.
    #       4.  Create a plot of frequency (%) vs time.
    #       5.  Depending on the frequency of each feature, different pharmacophore
    #           models can be derived.
    #           
    #           Features that have really low frequencies can be disregarded. (maybe < 5%)
    #           Features can be grouped by time-steps. Features with high frequency in some time-step
    #           can be used to build a pharmacophore model, while the ones that have low frequency in that
    #           time-step are not part of that pharmacophore model.           
    #
    #       6. How do we obtain the coordinates of the pharmacophore models? 
    #           Perhaps we could take the average for those timesteps.
    #
    #       7. This method can give one or more pharmacophores.
    #
    #       See: Wieder, Marcus, Ugo Perricone, Thomas Seidel, and Thierry Langer. "Pharmacophore models 
    #           derived from molecular dynamics simulations of protein-ligand complexes: A case study." 
    #           Natural product communications 11, no. 10 (2016): 1934578X1601101019.
    #
    # Common-Hits-Approach (CHA):
    #       Step 1. Getting Representative Pharmacophore Models.
    #
    #       1. Derive a pharmacophore model for each frame of the trajectory.
    #       2. Each pharmacophore can be respresented with a coordinates matrix and a feature vector, that is a vector with just
    #          the type of the features and no coordinates. 
    #       3. All pharmacophores with the same feature vector are reduced to only one. This is called the representative
    #          pharmacophore model or RPM.
    #       4. Get all unique pharmacophore features in all the trajectory (just feature type and atoms that are part of the feature).
    #          Build a vector with this features, unique features vector.
    #       5. Get all distint feature vectors by observing the feature vectors at each time step. 
    #          
    #          The feature vectors can be built as bit vectors with 1 corresponding to features
    #          that are present in the unique features vector and 0 to features which are not.
    #          
    #          Consider only vectors with an appereance count >= 2.
    #       6. For each time step a distinct feature vector appears, compute the conformational energy of the ligand. Save
    #          all of these energies and pick the median energy.
    #       7. Get the pharmacophore coordinates of the ligand with median energy.
    #       8. Build the RPM with the coordinates and the distinct feature vector
    #       9. Go back to 6 and repeat for all distinct feature vectors to get all RPMs for the system. 
    #
    #       Step 2. CHA
    #       1. Perform screening with each RPM.
    #       2. Compare lists of hits to find best phamracophore.
    # 
    # 
    # 
    # 
    # 
    #       See: Wieder, Marcus, Arthur Garon, Ugo Perricone, Stefan Boresch, Thomas Seidel, Anna Maria Almerico, 
    #            and Thierry Langer. "Common hits approach: combining pharmacophore modeling and molecular dynamics 
    #            simulations." Journal of chemical information and modeling 57, no. 2 (2017): 365-385       
    #