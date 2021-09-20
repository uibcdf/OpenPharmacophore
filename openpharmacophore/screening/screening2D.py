from rdkit.Chem.Pharm3D.Pharmacophore import Pharmacophore
from openpharmacophore.screening.screening import VirtualScreening
from rdkit import Chem, DataStructs
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D
from rdkit.Chem.Pharm2D.Generate import Gen2DFingerprint
import numpy as np
import pandas as pd
import json

class VirtualScreening2D(VirtualScreening):
    """ Class to perform virtual screening using pharmacophore fingerprints.

    Parameters
    ----------
        molecule: rdkit.Chem.mol
            The molecule whose pharmacophoric fingerprint will be comupted and
            used as query for screening.

        similarity: str (optional)
            Similarity measure that will be used to compare fingerprints.
            Defaults to tanimoto.
        
        sim_cutoff: float between 0 and 1
            Cutoff value from which a molecule is considered similar to the 
            query pharmacophore. Molecules below this value will not be kept

    Attributes
    ----------
        similarity_fn: str
            The similarity function that will be used to compare fingerprints.

    This class accepts a single molecule as query cause rdkit pharmacophore fingerprints 
    are limited to a single molecule.

    In the future this should be updated to create fingerprints for a consensus pharmacophore. 

    """
    def __init__(self, molecule, similarity="tanimoto", sim_cutoff=0.6):
        super().__init__(pharmacophore=self._get_pharmacophore_fingerprint(molecule))
        
        if similarity != "tanimoto" and similarity != "dice":
            raise NotImplementedError
        
        if sim_cutoff < 0 and sim_cutoff > 1:
            raise ValueError("Similarity cutoff value must lie between 0 and 1")

        self.n_similar_mols = self.n_matches
        self.scoring_metric = "Similarity"
        self.similarity_cutoff = sim_cutoff
        self.similiarity_fn = similarity
        self.similar_mols = self.matches

        self._factory = Gobbi_Pharm2D.factory
        self._screen_fn = self._fingerprint_similarity
    
    def print_report(self):
        """ Prints a summary report of the screening.
        """
        # TODO: Put this method in the base class to avoid repetition.
        report_str = "Virtual Screening Results\n"
        report_str += "-------------------------\n"
        report_str += "\nMolecules scanned: " 
        report_str += "{:,}".format(self.n_molecules).rjust(36)
        report_str += "\nMolecules matched to pharmacophore: " 
        report_str += str(self.n_similar_mols).rjust(19)
        report_str += "\nMolecules that didn't match the pharmacophore: " 
        report_str += "{:,}".format(self.n_fails).rjust(8)
        if self.n_similar_mols > 0:
            report_str += "\nLowest  similarity value: " + str(round(self.similar_mols[0][0], 4)).rjust(10)
            report_str += "\nHighest similarity value: " + str((round(self.similar_mols[-1][0], 4))).rjust(10)
            # Calculate mean similarity
            mean = 0
            N = self.n_similar_mols
            for i in range(N):
                mean += self.similar_mols[i][0]
            mean /= N
            report_str += "\nAverage similarity value: " + str(round(mean, 4)).rjust(10)         
            # Print top 5 molecules or less if there are less than 5
            if self.n_similar_mols < 5:
                n_top_mols = min(self.n_similar_mols, 5)
            else:
                n_top_mols = 5                
            report_str += "\n\nTop {} molecules:\n".format(n_top_mols)
            report_str += "\nZincID " + "Similarity".rjust(7)
            report_str += "\n-------  " + " ------\n"
            for i in range(n_top_mols):
                report_str += str(self.similar_mols[i][-1]) + "   "
                report_str += str(round(self.similar_mols[i][0], 4)) + "\n"
        print(report_str)

    def save_results_to_file(self, file_name):
        """Save the results of the screening to a file. The file contains the 
           matched molecules ids, smiles and similarity value. File can be 
           saved as csv or json.

           Parameters
           ----------
           file_name: str
                Name of the file

           Notes
           -----
           Does not return anything. A new file is written.
        """
        # TODO: Put this method in the base class to avoid repetition.
        file_format = file_name.split(".")[-1]

        similarity = [i[0] for i in self.similar_mols]
        smiles = [Chem.MolToSmiles(i[1]) for i in self.similar_mols]
        ids = [i[2] for i in self.similar_mols]

        results = {
            "mol_id": ids,
            "smiles": smiles,
            self.similiarity_fn: similarity,
        }

        if file_format == "csv":
            df = pd.DataFrame().from_dict(results)
            df.to_csv(file_name, index=False)
        elif file_format == "json":
            json_str = json.dumps(results)
            with open(file_name, "w") as f:
                f.write(json_str)
        else:
            raise NotImplementedError

    def _get_pharmacophore_fingerprint(self, molecule):
        """ Compute pharmacophore fingreprint for a single molecule

            Parameters
            ----------
            molecule: rdkit.Chem.mol
            
            Returns
            -------
            fingerprint: rdkit.DataStructs.SparseBitVect
        """
        factory = Gobbi_Pharm2D.factory
        fingerprint = Gen2DFingerprint(molecule, factory)
        return fingerprint
    
    def _fingerprint_similarity(self, molecules):
        """ Compute fingerprints and similarity values for a list
            of molecules. 

        Parameters
        ----------
        molecules: list of rdkit.Chem.mol
            List of molecules whose similarity to the pharmacophore 
            fingerprint will be calculated.
        
        Notes
        -------
        Does not return anything. Attributes are updated accordingly.

        """
        fps = [Gen2DFingerprint(mol, self._factory) for mol in molecules]
        self.n_molecules += len(fps)

        query_pharmacophore = self.pharmacophore
        if self.similiarity_fn == "tanimoto":
            sim_values = DataStructs.BulkTanimotoSimilarity(query_pharmacophore, fps)
        elif self.similiarity_fn == "dice":
            sim_values = DataStructs.BulkDiceSimilarity(query_pharmacophore, fps)
        else:
            raise NotImplementedError
        
        # TODO: try with a for loop to see if speed increases
        sim_values = np.array(sim_values)
        sim_mols_inx = (sim_values >= self.similarity_cutoff).nonzero()[0].tolist()

        similar_mols = [molecules[i] for i in sim_mols_inx]
        sim_values_filtered = sim_values[sim_mols_inx].tolist()
        mols_ids = [mol.GetProp("_Name") for mol in similar_mols]

        self.similar_mols += list(zip(sim_values_filtered, similar_mols, mols_ids))
        self.n_similar_mols += len(similar_mols)

        self.n_fails += self.n_molecules - self.n_similar_mols

