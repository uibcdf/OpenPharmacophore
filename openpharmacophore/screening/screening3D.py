from openpharmacophore.screening.screening import RetrospectiveScreening, VirtualScreening
from openpharmacophore.screening.alignment import apply_radii_to_bounds, transform_embeddings
from rdkit import RDConfig, Chem
from rdkit.Chem import ChemicalFeatures, rdDistGeom
from rdkit.Chem.Pharm3D import EmbedLib
import bisect
from operator import itemgetter
import os
    
class VirtualScreening3D(VirtualScreening):
    """ Class to perform virtual screening with a pharmacophore by 3D alignment of the molecules
        to the pharmacophore.

        Inherits from VirtualScreening class. 

    Code adapted from:
    https://github.com/rdkit/UGM_2016/blob/master/Notebooks/Stiefl_RDKitPh4FullPublication.ipynb

    Parameters
    ----------
    pharmacophore: openpharmacophore.Pharmacophore
        The pharmacophore that will be used to screen the database.

    Attributes
    ----------

    aligned_mols : list of 3-tuples (float, str, rdkit.Chem.mol)
        List of molecules that match the pharmacophore. Each tuple is formed by SSD value,
        the molecule id, and the molecule object.

    """

    def __init__(self, pharmacophore):
        super().__init__(pharmacophore)
        self.aligned_mols = self.matches 
        self.scoring_metric = "SSD"
        self._screen_fn = self._align_molecules
        
    def _align_molecules(self, molecules, verbose=0):
        """ Align a list of molecules to a given pharmacophore.

        Parameters
        ----------
        molecules: list of rdkit.Chem.mol
            List of molecules to align.

        verbose: int
            Level of verbosity
        
        Notes
        -------
        Does not return anything. The attributes aligned_mols, n_mathces,
        and n_fails are updated.

        """
        self.n_molecules += len(molecules)

        rdkit_pharmacophore, radii = self.pharmacophore.to_rdkit()
        apply_radii_to_bounds(radii, rdkit_pharmacophore)

        fdef = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
        featFactory = ChemicalFeatures.BuildFeatureFactory(fdef)

        for i, mol in enumerate(molecules):

            if verbose == 1 and i % 100 == 0 and i != 0:
                print(f"Screened {i} molecules. Number of matches: {self.n_matches}; Number of fails: {self.n_fails}")

            bounds_matrix = rdDistGeom.GetMoleculeBoundsMatrix(mol)
            # Check if the molecule features can match with the pharmacophore.
            # TODO: replace this function with a custom one that can take other feature definitions
            can_match, all_matches = EmbedLib.MatchPharmacophoreToMol(mol, featFactory, rdkit_pharmacophore)
            # all_matches is a list of tuples where each tuple contains the chemical features
            if can_match:
                # Match the molecule to the pharmacophore without aligning it
                failed, bounds_matrix_matched, matched_mols, match_details = EmbedLib.MatchPharmacophore(all_matches, 
                                                                                                bounds_matrix,
                                                                                                rdkit_pharmacophore, 
                                                                                                useDownsampling=True)
                if failed:
                    if verbose == 2:
                        print(f"Couldn't embed molecule {i}")
                    self.n_fails += 1
                    continue
            else:
                if verbose == 2:
                    print(f"Couldn't match molecule {i}")
                self.n_fails += 1
                continue
            atom_match = [list(x.GetAtomIds()) for x in matched_mols]
            try:
                mol_H = Chem.AddHs(mol)
                # Embed molecule onto the pharmacophore
                # embeddings is a list of molecules with a single conformer
                b_matrix, embeddings, num_fail = EmbedLib.EmbedPharmacophore(mol_H, atom_match, rdkit_pharmacophore, count=10)
            except Exception as e:
                if verbose == 2:
                    print(e)
                    print (f"Bounds smoothing failed for molecule {i}")
                self.n_fails += 1
                continue
            # Align embeddings to the pharmacophore 
            SSDs = transform_embeddings(rdkit_pharmacophore, embeddings, atom_match) 
            best_fit_index = min(enumerate(SSDs), key=itemgetter(1))[0]

            mol_id = mol.GetProp("_Name")
            matched_mol = (SSDs[best_fit_index], mol_id, embeddings[best_fit_index])
            # Append to list in ordered manner
            try:
                bisect.insort(self.aligned_mols, matched_mol) 
                self.n_matches += 1
            except:
                # Case when a molecule is repeated. It will throw an error since bisect
                # cannot compare molecules.
                self.n_molecules -= 1
                continue

class RetrospectiveScreening3D(RetrospectiveScreening):
    """ Class for performing retrospective virtual screening by 
        3D alignment of the molecules to the pharmacophore.

    Parameters
    ----------

    Attributes
    ----------

    """
    def __init__(self, pharmacophore):
        super().__init__(pharmacophore)
