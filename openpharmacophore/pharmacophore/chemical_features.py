# Openpharmacophore
from openpharmacophore._private_tools.exceptions import NoConformersError, OpenPharmacophoreTypeError
from openpharmacophore import PharmacophoricPoint
from openpharmacophore.utils.bisection import insort_right
# Third Party
import numpy as np
from rdkit import RDConfig, Chem
from rdkit.Chem import ChemicalFeatures
import pyunitwizard as puw
# Standard Library
import os
import pkg_resources
from  typing import Dict, Callable, Optional, List, Sequence, Tuple

rdkit_to_oph = {
        # To map rdkit feature names to openpharmacophore ones
        "Acceptor": "hb acceptor",
        "Donor": "hb donor",
        "Aromatic": "aromatic ring",
        "Hydrophobe": "hydrophobicity",
        "PosIonizable": "positive charge",
        "NegIonizable": "negative charge",
    }

def oph_featuredefinition() -> Dict[str, str]:
    """ Load default openpharmacophore feature definition.
    
        Returns
        -------
        dict
            Dictionary with chemical feature definitions (SMARTS strings).
    """
    feat_file = pkg_resources.resource_filename("openpharmacophore", "./data/smarts_features.txt")
    return load_smarts_fdef(feat_file)

def rdkit_featuredefinition() -> ChemicalFeatures.MolChemicalFeatureFactory:
    """ Loads rdkit chemical feature factory.
    
        Returns
        -------
        rdkit.Chem.rdMolChemicalFeatures.MolChemicalFeatureFactory
            The feature factory.

    """
    fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
    return ChemicalFeatures.BuildFeatureFactory(fdefName)

def load_smarts_fdef(file_name: str) -> Dict[str, str]:
    """ Load custom chemical feature definitions from a txt file.
    
        Feature definitions are SMART strings with their respective chemical feature
        name. 

        Parameters
        ----------
        file_name : str
            Name of the file containing the smarts feature definitions
        
        Returns
        -------
        features : dict
            Dictionary which keys are SMARTS strings and values are feature names

        Note
        -----
        The file must contain a SMARTS string follwed by the feature name. 
        
        Example:

        # Aromatic
        a1aaaaa1 Aromatic
        a1aaaa1 Aromatic

        Lines started with # are considered as comments
    """
    features = {} 
    # Load custom features file
    with open(file_name, "r") as file: 
        for line in file:
            if line[0] == "#" or line[0] == '\n':
                continue
            line = line.split(' ')
            feat_def = line[0] # The smarts string
            feat_name = line[1].rstrip()
            try:
                feat_name = rdkit_to_oph[feat_name]
            except KeyError:
                continue
            features[feat_def] = feat_name

    return features

class PharmacophoricPointExtractor():
    """ Class to extract pharmacophoric points from a ligand."""
    
    
    def __init__(self, featdef: Callable = oph_featuredefinition(), default_radius: float = 1.0, 
                directionality: bool = False, features: Optional[List[str]] = None) -> None:
        self.featdef = featdef
        self.default_radius = default_radius
        self.directionality = directionality
        if features is None:
            self.features = ['hb acceptor', 'aromatic ring', 'hb donor', 'hydrophobicity', 'positive charge', 'negative charge']
        else:
            self.features = features
    
    def extract_features(self, ligand: Chem.Mol, conformer_index: int) -> List[PharmacophoricPoint]:
        """ Extract the pharmacophoric points from a ligand.
        
            Parameters
            ----------
            ligand : rdkit.Chem.Mol
                A ligand
            
            conformer_index : int
                The conformer whose coordinates will be used to obtain the pharmacophoric
                points.
            
            Returns
            -------
            pharmacophoric_points : list of openpharmacophore.PharmacophoricPoint
                List of pharmacophoric points.
                
        """
        if ligand.GetNumConformers() == 0:
            raise NoConformersError
        
        if not isinstance(conformer_index, int):
            raise OpenPharmacophoreTypeError("conformer_index must be of type int")
        
        pharmacophoric_points = []
        
        if isinstance(self.featdef, dict):
            for smarts_pattern, feat_name in self.featdef.items():
                if feat_name not in self.features:
                    continue
                pattern = Chem.MolFromSmarts(smarts_pattern)  
                atom_indices = ligand.GetSubstructMatch(pattern)
                if len(atom_indices) == 0:
                    continue
                
                pharmacophoric_pnt = self.get_pharmacophoric_point(ligand, feat_name, atom_indices, 
                                                                   conformer_index, self.default_radius,
                                                                   self.directionality)
                insort_right(pharmacophoric_points, pharmacophoric_pnt, key=lambda p: p.short_name)
        else:
            # Use rdkit feature factory
            chemical_features = self.featdef.GetFeaturesForMol(ligand)
            for feature in chemical_features:
                feat_name = feature.GetFamily()
                feat_name = rdkit_to_oph[feat_name]
                if feat_name not in self.features:
                    continue
                atom_indices = feature.GetAtomIds()
                pharmacophoric_pnt = self.get_pharmacophoric_point(ligand, feat_name, atom_indices, 
                                                                   conformer_index, self.default_radius,
                                                                   self.directionality)
                insort_right(pharmacophoric_points, pharmacophoric_pnt, key=lambda p: p.short_name)
        
        return pharmacophoric_points
    
    @staticmethod
    def get_pharmacophoric_point(ligand: Chem.Mol, feat_name: str, atom_indices: Sequence, 
                                conformer_index: int, radius: float, directionality: bool) -> PharmacophoricPoint:
        """ Obtain the coordinates and if specified the direction vector and return a pharmacophoric point.
        
            Parameters
            ----------
            ligand : rdkit.Chem.Mol
                A ligand
            
            conformer_index : int
                The conformer whose coordinates will be used to obtain the pharmacophoric
                points.
            
            radius : float
                Lenght of the radius in angstroms of the parmacohporic point.
                
            directionality : bool
                Whether to compute the direction vectgor of that point.
            
            Returns
            -------
            openpharmacophore.PharmacophoricPoint
                A pharmacophoric point.
                
        """
        if len(atom_indices) > 1: 
            # Find the centroid
            # Aromatic, hydrophobic, positive or negative feature
            coords = PharmacophoricPointExtractor._feature_centroid(ligand, atom_indices, conformer_index)
            # Find direction vector
            if directionality:
                direction = PharmacophoricPointExtractor._aromatic_direction_vector(ligand, atom_indices, conformer_index)
            else:
                direction = None
        else:
            # Find the centroid
            # Donor or acceptor feature
            position = ligand.GetConformer(conformer_index).GetAtomPosition(atom_indices[0]) 
            coords = np.zeros((3,))
            coords[0] = position.x
            coords[1] = position.y
            coords[2] = position.z
            # Find direction vector
            if directionality:
                direction = PharmacophoricPointExtractor._donor_acceptor_direction_vector(ligand, atom_indices[0], coords, conformer_index)
            else:
                direction = None
                
        return PharmacophoricPoint(
                feat_type=feat_name,
                center=puw.quantity(coords, "angstroms"),
                radius=puw.quantity(radius, "angstroms"),
                direction=direction,
                atom_indices=atom_indices
            )
    
    @staticmethod                
    def _feature_centroid(molecule: Chem.Mol, atom_indxs: Tuple[int], conformer_index: int) -> np.ndarray:
        """
            Get the 3D coordinates of the centroid of a feature that encompasses more than 
            one atom. This could be aromatic, hydrophobic, negative and positive features

            Parameters
            ----------
            molecule : rdkit.Chem.Mol
                    Molecule that contains the feature which centroid will be computed

            atom_indxs : tuple of int
                    Indices of the atoms that belong to the feature

            conformer_index : int 
                    Index of the conformer for which the feature centroid will be computed

            Returns
            -------
            centroid : numpy.ndarray of shape (3, )
                Array with the coordinates of the centroid of the feature.

        """
        
        n_atoms = len(atom_indxs)
        coords = np.zeros((n_atoms, 3))
        for j, idx in enumerate(atom_indxs):
                position = molecule.GetConformer(conformer_index).GetAtomPosition(idx)
                coords[j, 0] = position.x
                coords[j, 1] = position.y
                coords[j, 2] = position.z
        
        centroid = coords.mean(axis=0)
    
        return centroid
    
    @staticmethod
    def _donor_acceptor_direction_vector(molecule: Chem.Mol, feat_type: str, atom_indx: int, 
                                        coords: np.ndarray, conformer_idx: int) -> np.ndarray:
        """
            Compute the direction vector for an H bond donor or H bond acceptor feature 

            Parameters
            ----------
            molecule : rdkit.Chem.rdchem.Mol
                    Molecule that contains the feature which direction vector will be computed.

            feat_type : str
                    Type of feature. Wheter is donor or acceptor.

            atom_indx : int
                    Index of the H bond acceptor or donor atom.

            coords : numpy.ndarray; shape(3,)
                    Coordiantes of the H bond acceptor or donor atom.

            conformer_idx : int 
                    Index of the conformer for which the direction vector will be computed.

            Returns
            -------
            direction : numpy.ndarray; shape(3,)
                    Coordinates of the direction vector.

        """
        direction = np.zeros((3,)) 
        atom = molecule.GetAtomWithIdx(atom_indx)
        for a in atom.GetNeighbors():
            if a.GetSymbol() == "H":
                continue
            position = molecule.GetConformer(conformer_idx).GetAtomPosition(a.GetIdx())
            direction[0] += position.x - coords[0]
            direction[1] += position.y - coords[1]
            direction[2] += position.z - coords[2]
        if feat_type == "Donor":
            direction = -direction
        return direction
        
    @staticmethod
    def _aromatic_direction_vector(molecule: Chem.Mol, atom_indxs: Tuple[int], conformer_idx: int) -> np.ndarray:
        """ Compute the direction vector for an aromatic feature. 

            Parameters
            ----------
            molecule : rdkit.Chem.rdchem.Mol
                    Molecule that contains the feature which direction vector will be computed.

            atom_indxs : tuple of int
                    Indices of the aromatic atoms.

            conformer_idx : int 
                    Index of the conformer for which the direction vector will be computed.

            Returns
            -------
            direction : numpy.ndarray; shape(3,)
                    Coordinates of the direction vector.

        """
        coords = np.zeros((3, 3)) # Take just the first three atoms
        for j, idx in enumerate(atom_indxs[0:3]):
                position = molecule.GetConformer(conformer_idx).GetAtomPosition(idx)
                coords[j, 0] = position.x
                coords[j, 1] = position.y
                coords[j, 2] = position.z
        
        # Find the vector normal to the plane defined by the three atoms
        u = coords[1, :] - coords[0, :]
        v = coords[2, :] - coords[0, :]
        direction = np.cross(u, v)

        return direction
    
    def __call__(self, ligand: Chem.Mol, 
                conformer_index: int) -> Callable[[Chem.Mol, int], List[PharmacophoricPoint]]:
        return self.extract_features(ligand, conformer_index)
    
    def __repr__(self) -> str:
         return f"{self.__class__.__name__}()"