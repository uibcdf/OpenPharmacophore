from rdkit import Chem

def load_mol2_file(fname):
    """Load molecules from a mol2 file
       
       code by https://chem-workflows.com/articles/2020/03/23/building-a-multi-molecule-mol2-reader-for-rdkit-v2/

    """
    molecules = []
    with open(fname, 'r') as f:
        doc = [line for line in f.readlines()]

    start = [index for (index,p) in enumerate(doc) if '@<TRIPOS>MOLECULE' in p]
    finish = [index-1 for (index,p) in enumerate(doc) if '@<TRIPOS>MOLECULE' in p]
    finish.append(len(doc))
    
    interval = list(zip(start, finish[1:]))
    
    for i in interval:
        block = ",".join(doc[i[0]:i[1]]).replace(',','')
        mol = Chem.MolFromMol2Block(block)
        molecules.append(mol)
    
    return molecules