from openpharmacophore._private_tools.exceptions import NoConformersError
from openpharmacophore.utils.centroid import feature_centroid
from openpharmacophore.utils.direction_vector import aromatic_direction_vector, donor_acceptor_direction_vector
from openpharmacophore.pharmacophoric_point import PharmacophoricPoint
import pyunitwizard as puw
from rdkit import RDConfig, Chem
from rdkit.Chem import ChemicalFeatures
import numpy as np
import os

def rdkit_points(ligands, radius, feat_list=None, direction_vector=False):
    """
        Get pharmacophoric points for a list of ligands using rdkit chemical feature definition. 

        Parameters
        ----------
        ligands: :obj: list of rdkit.Chem.rdchem.Mol
            List of ligands. Each ligand needs to have at least one conformer.
        
        radius: float
            Lenght of the radius of the parmacohporic points. Required if point type is 'spheres' or
            'spheres_vectors'.

        feat_list: list of str
            List of features that will be used to compute the pharmacophore.

        direction_vector: bool
            If true aromatic, donor and acceptors points will have direction.
        
        Returns
        -------
        points: dict
            Nested dictionary with the following structure:

            {ligand_id: 
                conformer_id: {
                    list of openpharmacophore.Pharmacophoric_Point
                    }
            }

            It stores pharmacophoric elements for each conformer of each ligand in the original
            ligand list.

    """
    fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    
    points = {}    
    for i, ligand in enumerate(ligands):

        n_conformers = ligand.GetNumConformers()
        if n_conformers == 0:
            raise NoConformersError(n_conformers)

        ligand_id = "ligand_" + str(i)
        points[ligand_id] = {}

        feats = factory.GetFeaturesForMol(ligand)
        
        for f in feats:
            
            feat_name = f.GetFamily()
            if feat_name not in feat_list:
                continue
            atom_idxs = f.GetAtomIds()
            
            for conformer_idx in range(n_conformers):
                if len(atom_idxs) > 1: 
                    # Find the centroid
                    coords = feature_centroid(ligand, atom_idxs, conformer_idx) # Aromatic, hydrophobic, positive or negative feature
                    # Find direction vector
                    if direction_vector:
                        direction = aromatic_direction_vector(ligand, atom_idxs, conformer_idx)
                    else:
                        direction = None
                else:
                    # Find the centroid
                    position = ligand.GetConformer(conformer_idx).GetAtomPosition(atom_idxs[0]) # Donor or acceptor feature
                    coords = np.zeros((3,))
                    coords[0] = position.x
                    coords[1] = position.y
                    coords[2] = position.z
                    # Find direction vector
                    if direction_vector:
                        direction = donor_acceptor_direction_vector(ligand, feat_name, atom_idxs[0], coords, conformer_idx)
                    else:
                        direction = None
                        
                point = rdkit_to_point(feat_name, coords, radius=radius, direction=direction, atom_indices=atom_idxs)
                conformer_id = "conformer_" + str(conformer_idx)

                if conformer_id not in points[ligand_id]:
                    points[ligand_id][conformer_id] = []
                points[ligand_id][conformer_id].append(point)
                
    return points

def custom_definition_points(ligands, radius, feat_list, feat_def, direction_vector=False):
    """
        Get pharmacophoric points for a list of ligands using custom smarts feature definition. 

        Parameters
        ----------
        ligands: :obj: list of rdkit.Chem.rdchem.Mol
            List of ligands. Each ligand needs to have at least one conformer.
        
        radius: float
            Lenght of the radius of the parmacohporic points. Required if point type is 'spheres' or
            'spheres_vectors'.

        feat_list: list of str
            List of features that will be used to compute the pharmacophore.
        
        feat_def: dict
            Definitions of the pharmacophoric points. 
            Dictionary which keys are SMARTS strings and values are feature names.
        
        direction_vector: bool
            If true aromatic, donor and acceptors points will have direction.

        Returns
        -------
        points: dict
            Nested dictionary with the following structure:

            {ligand_id: 
                conformer_id: {
                    list of openpharmacophore.Pharmacophoric_Point
                    }
            }
            
            It stores pharmacophoric elements for each conformer of each ligand in the original
            ligand list.

    """
    
    points = {}

    for i, ligand in enumerate(ligands):
        n_conformers = ligand.GetNumConformers()
        if n_conformers == 0:
            raise NoConformersError(n_conformers)
        
        ligand_id = "ligand_" + str(i)
        points[ligand_id] = {}

        for feat, feat_name in feat_def.items():
            if feat_name not in feat_list:
                continue
            pattern = Chem.MolFromSmarts(feat)  
            atom_idxs = ligand.GetSubstructMatch(pattern)
            if len(atom_idxs) == 0:
                continue
            for conformer_idx in range(n_conformers):
                if len(atom_idxs) > 1: 
                    # Find the centroid
                    coords = feature_centroid(ligand, atom_idxs, conformer_idx) # Aromatic, hydrophobic, positive or negative feature
                    # Find direction vector
                    if direction_vector:
                        direction = aromatic_direction_vector(ligand, atom_idxs, conformer_idx)
                    else:
                        direction = None
                else:
                    # Find the centroid
                    position = ligand.GetConformer(conformer_idx).GetAtomPosition(atom_idxs[0]) # Donor or acceptor feature
                    coords = np.zeros((3,))
                    coords[0] = position.x
                    coords[1] = position.y
                    coords[2] = position.z
                    # Find direction vector
                    if direction_vector:
                        direction = donor_acceptor_direction_vector(ligand, atom_idxs[0], coords, conformer_idx)
                    else:
                        direction = None

                point = rdkit_to_point(feat_name, coords, radius=radius, direction=direction, atom_indices=atom_idxs)
                conformer_id = "conformer_" + str(conformer_idx)

                if conformer_id not in points[ligand_id]:
                    points[ligand_id][conformer_id] = []
                points[ligand_id][conformer_id].append(point)

    return points


def ligands_pharmacophoric_points(ligands, radius, feat_list=None, feat_def=None):
    """
        Get pharmacophoric points for each ligand in a list of ligands. If a ligand has 
        more than one conformer, pharamcophoric points will be computed for each one.  

        Parameters
        ----------
        ligands: list of rdkit.Chem.Mol or a single rdkit.Chem.Mol 
            List of ligands or a single ligand. Each ligand needs to have at least one conformer.
        
        radius: float
            Lenght of the radius of the parmacohporic points. Required if point type is 'spheres' or
            'spheres_vectors'.

        feat_list: list of str (optional)
            List of features that will be used to compute the pharmacophore.

        feat_def: dict (optional)
            Definitions of the pharmacophoric points. 
            Dictionary which keys are SMARTS strings and values are feature names.

        Returns
        -------
        points: dict
            Nested dictionary with the following structure:

            {ligand_id: 
                conformer_id: {
                    list of openpharmacophore.Pharmacophoric_Point
                    }
            }
            
            It stores pharmacophoric elements for each conformer of each ligand in the original
            ligand list.

    """
    if isinstance(ligands, Chem.rdchem.Mol): # Check if it's a single molecule
        ligands = [ligands]
    elif not isinstance(ligands, list):
        ligands = list(ligands)
    
    if not feat_list:
        feat_list = ['Acceptor', 'Aromatic', 'Donor', 'Hydrophobe', 'PosIonizable', 'NegIonizable']
    
    if feat_def is None: # If no feature definition is given use rdkit one
        points = rdkit_points(ligands=ligands, radius=radius, feat_list=feat_list)
    else:
        points = custom_definition_points(ligands=ligands, radius=radius, feat_list=feat_list, feat_def=feat_def)
    
    return points

    
def rdkit_to_point(feat_name, coords, radius=None, direction=None, atom_indices=None):
    """Transform an rdkit feature point to an openpharmacophore pharmacophoric point.

        Parameters
        ----------
        feat_name: str
            rdkit name of the feature point.

        coords: numpy.ndarray; shape: (3, )
            3D coordinates of the centroid of the feature.
        
        radius: float
            Lenght of the radius of the parmacohporic points. 

        direction: list, tuple, numpy.ndarray; shape:(3,)
            Unit vector. 

        Returns
        -------
        point: an openpharmacophore.pharmacophoric_point.PharmacophoricPoint
    
    """

    points = {
        "Acceptor": "hb acceptor",
        "Donor": "hb donor",
        "Aromatic": "aromatic ring",
        "Hydrophobe": "hydrophobicity",
        "PosIonizable": "positive charge",
        "NegIonizable": "negative charge",
    }

    point = PharmacophoricPoint(
        feat_type=points[feat_name],
        center=puw.quantity(coords, "angstroms"),
        radius=puw.quantity(radius, "angstroms"),
        direction=direction,
        atoms_inxs=atom_indices
    )

    return point