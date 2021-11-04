from openpharmacophore.utils.align_ligands import align_set_of_ligands
from openpharmacophore.utils.centroid import feature_centroid
from openpharmacophore.utils.ligand_features import rdkit_to_point
import numpy as np
from sklearn.cluster import DBSCAN
from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalFeatures
import os

def get_feature_clusters(feat_coords, eps, min_samples):
    """
    Get clusters of features for a ligand based pharmacophore.

    Parameters
    ----------

    feat_coords : dict
        Dictionary containing 3D coordinates for each feature type. Dictionary keys 
        are feature name and values are an numpy array of coordinates 
        
    eps : float
        The maximum distance between two pharmacophoric points for one to be considered 
        as in the neighborhood of the other. (Default: 2)

    min_samples : float
        Percentages of ligands that must contain a pharmacophoric point to be considered as a core point. 
        Must be a number between 0 and 1
    
    Returns
    ----------

    clusters : dict
        Dictionary with centroid of each cluster of features. Keys are feature name
        and values is a list of coordinates

    """
    
    clusters = {}
    for feat, coords in feat_coords.items():
        db_scan = DBSCAN(eps=eps, min_samples=min_samples).fit(coords)
        
        labels = db_scan.labels_
        core_samples_mask = np.zeros_like(labels, dtype=bool)
        core_samples_mask[db_scan.core_sample_indices_] = True
        
        centroids = []
        unique_labels = set(labels)
        for k in unique_labels:
            if k == -1:
                continue
            class_member_mask = (labels == k)
            cluster = feat_coords[feat][class_member_mask & core_samples_mask]
            cluster_centroid = cluster.mean(axis=0)
            centroids.append(cluster_centroid)
        
        clusters[feat] = centroids
        
    return clusters

def dbscan_pharmacophore(ligands, radius=1, eps=2, min_samples=0.75, feat_list=None, feat_def=None):
    """
    Compute a ligand based pharmacophore from a list of ligands, using a density based 
    clustering algorithm.
    
    Parameters
    ----------

    ligands : list of rdkit.Chem.rdchem.Mol or rdkit.Chem.SmilesMolSupplier or rdkit.Chem.SDMolSupplier
            List of ligands.

    radius : float, optional
        Lenght of the radius in angstroms of the parmacohporic points (Default: 1)

    eps : float, optional
        The maximum distance between two pharmacophoric points for one to be considered 
        as in the neighborhood of the other (default=2).

    min_samples : float
        Percentages of ligands that must contain a pharmacophoric point to be considered as a core point. 
        Must be anumber betwwen 0 and 1 (default=0.75).
    
    feat_list : list of str, optional
            List of features that will be used to compute the pharmacophore.
            
    feat_def : dict
            Definitions of the pharmacophoric points. 
            Dictionary which keys are SMARTS strings and values are feature names.

    Returns
    --------

    pharmacophoric_points : list of openpharmacophore.PharmacophoricPoints
        The pharmacophoric points of the common pharmacophore.

    aligned_ligands: list of rdkit.Chem.Mol
        A list containing the aligned ligands.

    """
    if min_samples < 0 or min_samples > 1:
        raise ValueError("min_samples must be a value between 0 and 1")
    
    aligned_ligands, _ = align_set_of_ligands(ligands)

    if feat_def is None: # If no feature definition is given use rdkit one
        fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
        factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    else:
        reverse_feat_dict = {}
        for smarts, feat_name in feat_def.items():
            if feat_name not in reverse_feat_dict:
                reverse_feat_dict[feat_name] = []
            reverse_feat_dict[feat_name].append(smarts)
    
    if not feat_list:
        feat_list = ['Acceptor', 'Aromatic', 'Donor', 'Hydrophobe', 'PosIonizable', 'NegIonizable']
    
    feat_coords = {}
    for feature in feat_list:
        feat_coords[feature] = []

        for ligand in aligned_ligands:
            if feat_def is None:
                feats = factory.GetFeaturesForMol(ligand, includeOnly=feature)
                for f in feats:
                    atom_idxs = f.GetAtomIds()
                    if len(atom_idxs) > 1: # Aromatic, hydrophobic, positive or negative feature
                        coords = feature_centroid(ligand, atom_idxs, 0)
                    else: # Donor or acceptor feature
                        position = ligand.GetConformer(0).GetAtomPosition(atom_idxs[0])
                        coords = np.zeros((3,))
                        coords[0] = position.x
                        coords[1] = position.y
                        coords[2] = position.z
                feat_coords[feature].append((coords.tolist()))
            else:
                for f in reverse_feat_dict[feature]:
                    pattern = Chem.MolFromSmarts(f)  
                    atom_idxs = ligand.GetSubstructMatch(pattern)
                    if len(atom_idxs) == 0:
                        continue
                    elif len(atom_idxs) == 1: # Donor or acceptor feature
                        position = ligand.GetConformer(0).GetAtomPosition(atom_idxs[0])
                        coords = np.zeros((3,))
                        coords[0] = position.x
                        coords[1] = position.y
                        coords[2] = position.z
                    elif len(atom_idxs) > 1: # Aromatic, hydrophobic, positive or negative feature
                        coords = feature_centroid(ligand, atom_idxs, 0)      
                 
                    feat_coords[feature].append((coords.tolist()))
        feat_coords[feature] = np.array(feat_coords[feature])

    feat_coords = {feature: coords for feature, coords in feat_coords.items() if len(coords) > 0} # remove features with no coordinates
    
    min_samples = int(min_samples * len(ligands)) 
    feature_clusters = get_feature_clusters(feat_coords, eps=eps, min_samples=min_samples)
   
    pharmacophoric_points = []
    for feature_type, coords in feature_clusters.items():
        for center in coords:            
            point = rdkit_to_point(feature_type, center, radius=radius, direction=None)
            pharmacophoric_points.append(point)

    return pharmacophoric_points, aligned_ligands