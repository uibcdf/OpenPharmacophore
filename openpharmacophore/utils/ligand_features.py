from rdkit import RDConfig, Chem
from rdkit.Chem import ChemicalFeatures
from openpharmacophore._private_tools.exceptions import NoConformersError, PointTypeError
from openpharmacophore.utils.centroid import feature_centroid
from openpharmacophore.utils.rdkit_to_point import rdkit_to_point
import numpy as np
import os

def ligands_pharmacophoric_points(ligands, radius, feat_list=None, feat_def='rdkit', 
                                point_type="spheres"):

    """
        Get pharmacophoric points for each ligand in a list of ligands. If a ligand has 
        more than one conformer, pharamcophoric points will be computed for each one.  

        Parameters
        ----------
        ligands: :obj: rdkit.Chem.rdmolfiles.SmilesMolSupplier or list of rdkit.Chem.rdchem.Mol
            List of ligands. Each ligand needs to have at least one conformer
        
        radius: float
            Lenght of the radius of the parmacohporic points. Required if point type is 'spheres' or
            'spheres_vectors'

        feat_list: list of str (optional)
            List of features that will be used to compute the pharmacophore

        feat_def: str
            Defenitions of the pharmacophoric points (Default: 'rdkit')
        
        point_type: str
            Type of pharmacophoric points to be returned

        Returns
        -------
        points: a list of openpharmacophore.pharmacophoric_elements

    """
    # TODO: add other feature definitions

    point_type_list = ["spheres", "spheres_vectors", "gaussian", "shapelet"]
    if point_type not in point_type_list:
        raise PointTypeError(point_type)

    if (point_type == "spheres" or point_type == "spheres_vectors") and radius is None:
        raise ValueError("Radius cannot be null if point type is spheres or spheres_vectors") 

    if feat_def == 'rdkit':
        fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
        factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    else:
        raise NotImplementedError

    if isinstance(ligands, Chem.rdchem.Mol): # Check if it's a single molecule
        ligands = [ligands]
    elif not isinstance(ligands, list):
        ligands = list(ligands)
    
    if not feat_list:
        feat_list = ['Acceptor', 'Aromatic', 'Donor', 'Hydrophobe', 'PosIonizable', 'NegIonizable']
    
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
                                    
                if len(atom_idxs) > 1: # Get the centroid of that feature
                    coords = feature_centroid(ligand, atom_idxs, conformer_idx)
                else:
                    position = ligand.GetConformer(conformer_idx).GetAtomPosition(atom_idxs[0])
                    coords = np.zeros((3,))
                    coords[0] = position.x
                    coords[1] = position.y
                    coords[2] = position.z

                # TODO: find direction vector for vectorial pharmacophoric points

                point = rdkit_to_point(feat_name, coords, radius=radius, point_type=point_type)

                conformer_id = "conformer_" + str(conformer_idx)

                if conformer_id in points[ligand_id]:
                    points[ligand_id][conformer_id].append(point)
                else:
                    points[ligand_id][conformer_id] = []

    return points

