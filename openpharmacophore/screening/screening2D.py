from openpharmacophore.screening.screening import VirtualScreening
from rdkit import DataStructs
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D
from rdkit.Chem.Pharm2D.Generate import Gen2DFingerprint
import bisect

class VirtualScreening2D(VirtualScreening):
    """ Class to perform virtual screening using pharmacophore fingerprints.
        Inherits from VirtualScreening class.

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
            query pharmacophore. Molecules below this value will not be kept.

    Attributes
    ----------

        similar_mols : list of 3-tuples (float, str, rdkit.Chem.mol)
            List of molecules that are considered similar enough to the
            phamracophore fingerprint.

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

        self.scoring_metric = "Similarity"
        self.similarity_cutoff = sim_cutoff
        self.similiarity_fn = similarity
        self.similar_mols = self.matches

        self._factory = Gobbi_Pharm2D.factory
        self._screen_fn = self._fingerprint_similarity
    
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
        for mol in molecules:
            self.n_molecules += 1
            fingerprint = Gen2DFingerprint(mol, self._factory)
            if self.similiarity_fn == "tanimoto":
                similarity = DataStructs.TanimotoSimilarity(self.pharmacophore, fingerprint)
            elif self.similiarity_fn == "dice":
                similarity = DataStructs.DiceSimilarity(self.pharmacophore, fingerprint)
            else:
                raise NotImplementedError
            
            if similarity >= self.similarity_cutoff:
                try:
                    mol_id = mol.GetProp("_Name")
                except:
                    mol_id = None
                matched_mol = (similarity, mol_id, mol)
                # Append to list in ordered manner
                try:
                    bisect.insort(self.similar_mols, matched_mol)
                    self.n_matches += 1
                except:
                    # Case when a molecule is repeated. It will throw an error since bisect
                    # cannot compare molecules.
                    self.n_molecules -= 1
                    continue
            else:
                self.n_fails += 1 


