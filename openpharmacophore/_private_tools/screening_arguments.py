from openpharmacophore import Pharmacophore, StructuredBasedPharmacophore, LigandBasedPharmacophore
from openpharmacophore._private_tools.exceptions import InvalidSimilarityFunction, OpenPharmacophoreValueError
from rdkit import DataStructs

def is_3d_pharmacophore(pharmacophore):
    """ Check whether a pharmacophore object is of type Phamracophore, StructuredBasedPharmacophore
        or LigandBasedPharmacophore
        
        Parameters
        ----------
        pharmacophore : obj
            The pharmacophore object.
    
        Returns
        -------
        bool
    """
    if (isinstance(pharmacophore, Pharmacophore) 
        or isinstance(pharmacophore, StructuredBasedPharmacophore)
        or isinstance(pharmacophore, LigandBasedPharmacophore)):
           return True
    return False

def check_virtual_screening_kwargs(**kwargs):
    """ Check if kwargs were passed to an instance of the virtual screening classes.
    
        Returns
        -------
        similiarity_fn : {DataStructs.TanimotoSimilarity, DataStructs.DiceSimilarity}
            The fingerprint similarity function.
        
        similarity_cutoff : float
            The similarity value form which a molecule is considered active.
    """
    if kwargs:
            if "similarity" in kwargs:
                sim = kwargs["similarity"]
                if sim == "tanimoto":
                    similiarity_fn = DataStructs.TanimotoSimilarity
                elif sim == "dice":
                    similiarity_fn = DataStructs.DiceSimilarity
                else:
                    raise InvalidSimilarityFunction(f"{sim} is not a valid similarity function")
                
            else:
                similiarity_fn = DataStructs.TanimotoSimilarity

            if "sim_cutoff" in kwargs:
                if kwargs["sim_cutoff"] < 0 and kwargs["sim_cutoff"] > 1:
                    raise OpenPharmacophoreValueError("Similarity cutoff value must lie between 0 and 1")
                similarity_cutoff = kwargs["sim_cutoff"]
            else:
                similarity_cutoff = 0.2
    else:
        similiarity_fn = DataStructs.TanimotoSimilarity
        similarity_cutoff = 0.2

    return similiarity_fn, similarity_cutoff
