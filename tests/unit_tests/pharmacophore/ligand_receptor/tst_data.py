# Add data used in tests such as molecules, trajectories, ...
from rdkit import Chem


estradiol_pdb_str = """
REMARK   1 CREATED WITH MDTraj 1.9.7, 2022-10-17
CRYST1  105.500  105.500  136.080  90.00  90.00 120.00 P 1           1 
MODEL        0
ATOM   5944  C1  EST A 600     104.106  17.203  24.775  1.00  0.00           C  
ATOM   5945  C2  EST A 600     102.995  17.834  25.370  1.00  0.00           C  
ATOM   5946  C3  EST A 600     101.695  17.355  25.120  1.00  0.00           C  
ATOM   5947  O3  EST A 600     100.598  17.990  25.704  1.00  0.00           O  
ATOM   5948  C4  EST A 600     101.506  16.240  24.274  1.00  0.00           C  
ATOM   5949  C5  EST A 600     102.621  15.588  23.660  1.00  0.00           C  
ATOM   5950  C6  EST A 600     102.371  14.379  22.735  1.00  0.00           C  
ATOM   5951  C7  EST A 600     103.644  13.753  22.086  1.00  0.00           C  
ATOM   5952  C8  EST A 600     104.898  13.873  22.953  1.00  0.00           C  
ATOM   5953  C9  EST A 600     105.178  15.388  23.261  1.00  0.00           C  
ATOM   5954  C10 EST A 600     103.957  16.078  23.918  1.00  0.00           C  
ATOM   5955  C11 EST A 600     106.462  15.459  24.125  1.00  0.00           C  
ATOM   5956  C12 EST A 600     107.711  14.803  23.508  1.00  0.00           C  
ATOM   5957  C13 EST A 600     107.463  13.343  23.124  1.00  0.00           C  
ATOM   5958  C14 EST A 600     106.170  13.270  22.242  1.00  0.00           C  
ATOM   5959  C15 EST A 600     106.228  11.821  21.792  1.00  0.00           C  
ATOM   5960  C16 EST A 600     107.701  11.713  21.263  1.00  0.00           C  
ATOM   5961  C17 EST A 600     108.494  12.719  22.135  1.00  0.00           C  
ATOM   5962  O17 EST A 600     109.610  12.027  22.746  1.00  0.00           O  
ATOM   5963  C18 EST A 600     107.379  12.449  24.419  1.00  0.00           C  
TER    5964      EST A 600
ENDMDL
CONECT    1    2   11
CONECT    2    1    3
CONECT    3    2    4    5
CONECT    4    3
CONECT    5    3    6
CONECT    6    5    7   11
CONECT    7    6    8
CONECT    8    7    9
CONECT    9    8   10   15
CONECT   10    9   11   12
CONECT   11    1    6   10
CONECT   12   10   13
CONECT   13   12   14
CONECT   14   13   15   18   20
CONECT   15    9   14   16
CONECT   16   15   17
CONECT   17   16   18
CONECT   18   14   17   19
CONECT   19   18
CONECT   20   14
END
    """

estradiol = Chem.MolFromPDBBlock(estradiol_pdb_str)
assert estradiol.GetNumAtoms() == 20


threonine_pdb_str = """
REMARK   1 CREATED WITH MDTraj 1.9.7, 2022-10-24
CRYST1   84.599   89.560   94.156  90.00  90.00  90.00 P 1           1 
MODEL        0
ATOM      1  N   THR A   2      44.235  80.308  18.419  1.00  0.00           N  
ATOM      2  CA  THR A   2      43.549  79.243  17.706  1.00  0.00           C  
ATOM      3  C   THR A   2      44.528  78.252  17.077  1.00  0.00           C  
ATOM      4  O   THR A   2      45.699  78.559  16.853  1.00  0.00           O  
ATOM      5  CB  THR A   2      42.608  79.792  16.611  1.00  0.00           C  
ATOM      6  OG1 THR A   2      43.375  80.468  15.608  1.00  0.00           O  
ATOM      7  CG2 THR A   2      41.586  80.762  17.208  1.00  0.00           C  
ENDMDL
END
"""

threonine = Chem.MolFromPDBBlock(threonine_pdb_str)
assert threonine.GetNumAtoms() == 7

del estradiol_pdb_str
del threonine_pdb_str
