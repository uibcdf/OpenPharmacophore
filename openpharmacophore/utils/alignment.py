from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdMolAlign
from openpharmacophore.utils.conformers import generate_conformers
import numpy as np
import copy

def align_set_of_ligands(ligands):
    """
        Align a set of ligands to each other
    """
    
    molecules = copy.deepcopy(list(ligands))
    molecules = [generate_conformers(mol, 100) for mol in molecules]

    crippen_contribs = [rdMolDescriptors._CalcCrippenContribs(mol) for mol in molecules]
    crippen_ref_contrib = crippen_contribs[0]
    crippen_prob_contribs = crippen_contribs

    ref_mol = molecules[0]
    probe_mols = molecules

    crippen_score = []
    aligned_molecules = []
    for idx, mol in enumerate(probe_mols):
        tempscore = []
        
        for cid in range(100):
            crippenO3A = rdMolAlign.GetCrippenO3A(mol, ref_mol, crippen_prob_contribs[idx], crippen_ref_contrib, cid, 0)
            crippenO3A.Align()
            tempscore.append(crippenO3A.Score())
            
        best = np.argmax(tempscore)
        mol_string = Chem.MolToMolBlock(mol, confId=int(best))
        temp_mol = Chem.MolFromMolBlock(mol_string, removeHs=False)
        
        crippen_score.append(tempscore[best])
        aligned_molecules.append(temp_mol)
    
    return aligned_molecules, crippen_score
        
    