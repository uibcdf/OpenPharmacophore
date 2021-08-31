from openpharmacophore.io.mol2 import load_mol2_file

from rdkit import RDConfig, Chem, Geometry
from rdkit.Chem import ChemicalFeatures, rdDistGeom, rdMolTransforms
from rdkit.Chem.Pharm3D import Pharmacophore, EmbedLib
from rdkit.Numerics import rdAlignment
import pyunitwizard as puw
from tqdm.auto import tqdm

from operator import itemgetter
import os
    
class VirtualScreening():
    """ Class to perform virtual screening with a pharmacophore

    Code adapted from:
    https://github.com/rdkit/UGM_2016/blob/master/Notebooks/Stiefl_RDKitPh4FullPublication.ipynb

    Parameters
    ----------

    Attributes
    ----------

    success : list of 2-tuples (float, rdkit.Chem.mol)
        List of molecules with its SSD value retrieved from the screen that match the pharmacophore

    n_success : int
        Number of molecules that match with the pharmacophore

    n_fails : int
         Number of molecules that do not match with the pharmacophore
    """

    def __init__(self):
        self.success = [] 
        self.n_success = 0 
        self.n_fails = 0
    
    def screen_db(self, pharmacophore, db, verbose=1):

        # TODO: 
        # Implement screening for multiple files
        # Fetch a database
        # Vectorial pharmacophoric points

        points = []
        radii = []

        rdkit_element_name = { # dictionary to map openpharmacophore feature names to rdkit feature names
        "aromatic ring": "Aromatic",
        "hydrophobicity": "Hydrophobe",
        "hb acceptor": "Acceptor",
        "hb donor": "Donor",
        "positive charge": "PosIonizable",
        "negative charge": "NegIonizable",
        }

        for element in pharmacophore.elements:
            feat_name = rdkit_element_name[element.feature_name]
            center = puw.get_value(element.center, to_unit="angstroms")
            center = Geometry.Point3D(center[0], center[1], center[2])
            points.append(ChemicalFeatures.FreeChemicalFeature(
                feat_name,
                center
            ))
            radius = puw.get_value(element.radius, to_unit="angstroms")
            radii.append(radius)

        rdkit_pharmacophore = Pharmacophore.Pharmacophore(points)
        apply_radii_to_bounds(radii, rdkit_pharmacophore)

        fdef = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
        featFactory = ChemicalFeatures.BuildFeatureFactory(fdef)

        molecules = load_db(db)
        self.success = []
        for i, mol in enumerate(tqdm(molecules)):

            if verbose == 1 and i % 100 == 0 and i != 0:
                print(f"Screened {i} molecules. Number of matches: {self.n_success}; Number of fails: {self.n_fails}")

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
            self.success.append((SSDs[best_fit_index], embeddings[best_fit_index]))
            self.n_success += 1

        
    def enrichment_plot():
        pass

    def enrichment_factor():
        pass

    def ROC_plot():
        pass

    def AUC():
        pass


def get_transform_matrix(align_ref, conf_embed, atom_match):

    """Get the transformation matrix for a conformer.

        Parameters
        ----------
        align_ref: list of rdkit.Geometry.Point3D
            list of pharmacophore reference points for the alignment.

        conf_embed: rdkit.Chem.Conformer
            The conformer embedding.

        atom_match: list of list
            List of list of atoms ids that match the pharmacophore.

        Returns
        -------
        a 2-tuple

        ssd: float
            SSD value for the alignment.

        transform_matrix: numpy.ndarray; shape(4, 4)
            The transform matrix.

        """
  
    align_probe = [] # probe points to align to reference points
    for match_ids in atom_match:
        dummy_point = Geometry.Point3D(0.0,0.0,0.0)
        for id in match_ids:
            dummy_point += conf_embed.GetAtomPosition(id)
        dummy_point /= len(match_ids)
        align_probe.append(dummy_point)
    
    ssd, transform_matrix = rdAlignment.GetAlignmentTransform(align_ref, align_probe)
    return ssd, transform_matrix

def transform_embeddings(pharmacophore, embeddings, atom_match):

    """Transform embeddings. Performs the alignment of the molecules 
        to the pharmacophore.

        Parameters
        ----------
        pharmacophore: rdkit.Chem.Pharm3D.Pharmacophore
            A pharmacophore object.

        embeddings: list of rdkit.Chem.Mol
            List of molecules with a single conformer.

        atom_match: list of list
            List of list of atoms ids that match the pharmacophore.

        Returns
        -------
        SSDs: list of float
            List of sum of square deviations (SSD) values for the alignments.

        """

    align_ref = [f.GetPos() for f in pharmacophore.getFeatures()]
    ssds = []
    for embedding in embeddings:
        conformer = embedding.GetConformer()
        ssd, transform_matrix = get_transform_matrix(align_ref, conformer, atom_match)
        # Transform the coordinates of the conformer
        rdMolTransforms.TransformConformer(conformer, transform_matrix)
        ssds.append(ssd)
    
    return ssds

def apply_radii_to_bounds(radii, pharmacophore):
    """
        Apply the radius of the pharmacophoric points to the bound
        matrix of the pharmacophore

        Parameters
        ----------
        radii: list of float
            list with the radius of each of the pharmacophoric points
        
        pharmacophore: rdkit.Chem.Pharm3D.Pharmacophore
            The pharmacophore wich bounds matriz will be modified
        
        Returns
        -------
        Nothing is returned. The pharmacophore bounds matrix is updated

    """
    for i in range(len(radii)):
        for j in range(i + 1, len(radii)):
            sum_radii = radii[i] + radii[j]
            pharmacophore.setLowerBound(i, j, max(pharmacophore.getLowerBound(i, j) - sum_radii, 0))
            pharmacophore.setUpperBound(i, j, pharmacophore.getUpperBound(i, j) + sum_radii)


def load_db(db):
    """
        Load a file of molecules of any format
        and return a list of rdkit molecules
    """
    if isinstance(db, str):

        fextension = db.split(".")[-1]
        
        if fextension == "smi":
            ligands = Chem.SmilesMolSupplier(db, delimiter=' ', titleLine=True)
        elif fextension == "mol2":
            ligands = load_mol2_file(db)
        elif fextension == "sdf":
            ligands = Chem.SDMolSupplier(db)
        else:
            raise NotImplementedError
        
        if len(ligands) == 0:
            raise Exception("Molecules couldn´t be loaded")
        
        return list(ligands)
    
    return db