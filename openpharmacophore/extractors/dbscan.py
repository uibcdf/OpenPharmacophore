from openpharmacophore.utils.alignment import align_set_of_ligands
from openpharmacophore.utils.rdkit_to_point import rdkit_to_point
import numpy as np
import os
from sklearn.cluster import DBSCAN
from rdkit import RDConfig
from rdkit.Chem import ChemicalFeatures
from openpharmacophore.utils.centroid import feature_centroid

def get_feature_clusters(feat_coords, eps, min_samples):
    """
    Get clusters of features for a ligand based pharmacophore
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

def dbscan_pharmacophore(ligands, radius=1, eps=2, min_samples=0.75, feat_list=None):
    """
    Compute a ligand based pharmacophore from a list of ligands, using a density based 
    clustering algorithm.
    """
    
    aligned_ligands, _ = align_set_of_ligands(ligands)

    fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    
    if not feat_list:
        feat_list = ['Acceptor', 'Aromatic', 'Donor', 'Hydrophobe', 'PosIonizable', 'NegIonizable']
    
    feat_coords = {}
        
    for feature in feat_list:
        feat_coords[feature] = []

        for ligand in aligned_ligands:
            feats = factory.GetFeaturesForMol(ligand, includeOnly=feature) 

            for f in feats:
                atom_idxs = f.GetAtomIds()
                if len(atom_idxs) > 1: # Get the centroid of that feature
                    coords = feature_centroid(ligand, atom_idxs, 0)
                else:
                    position = ligand.GetConformer(0).GetAtomPosition(atom_idxs[0])
                    coords = np.zeros((3,))
                    coords[0] = position.x
                    coords[1] = position.y
                    coords[2] = position.z

                feat_coords[feature].append(coords.tolist())
            
        feat_coords[feature] = np.array(feat_coords[feature])
        
    min_samples = int(min_samples * len(ligands)) 
    feature_clusters = get_feature_clusters(feat_coords, eps=eps, min_samples=min_samples)

    pharmacophoric_points = []
    for feature_type, coords in feature_clusters.items():
        for center in coords:
            point = rdkit_to_point(feature_type, center, radius=radius)
            pharmacophoric_points.append(point)

    return pharmacophoric_points, aligned_ligands