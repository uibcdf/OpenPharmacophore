import numpy as np
import pyunitwizard as puw

from openpharmacophore.molecular_systems.convert import mol_to_topology
from openpharmacophore.io.pdb import write_pdb_block


def test_save_topology_to_pdb(dialanine):
    topology = mol_to_topology(dialanine)
    coords = puw.quantity(np.array([[
        [46.411, 13.580, 27.508],
        [46.204, 13.395, 26.090],
        [47.270, 14.190, 25.332],
        [47.947, 13.654, 24.448],
        [44.817, 13.873, 25.699],
        [50.424, 13.298, 22.736],
        [50.417, 13.697, 21.333],
        [51.805, 14.263, 21.035],
        [52.408, 13.939, 20.017],
        [49.360, 14.751, 21.095],
        [45.750, 12.989, 28.032],
        [47.376, 13.314, 27.752],
        [46.286, 12.328, 25.836],
        [47.435, 15.250, 25.577],
        [44.230, 14.075, 26.607],
        [44.901, 14.794, 25.104],
        [44.316, 13.096, 25.104],
        [50.409, 14.138, 23.332],
        [49.593, 12.722, 22.933],
        [50.187, 12.844, 20.678],
        [52.268, 14.968, 21.741],
        [48.988, 15.122, 22.061],
        [49.796, 15.585, 20.525],
        [48.527, 14.313, 20.525],
    ]]), "angstroms")
    pdb_str = write_pdb_block(topology, coords)

    expected = """MODEL        0
ATOM      1  N   ALA A   0      46.411  13.580  27.508  1.00  0.00           N  
ATOM      2  CA  ALA A   0      46.204  13.395  26.090  1.00  0.00           C  
ATOM      3  C   ALA A   0      47.270  14.190  25.332  1.00  0.00           C  
ATOM      4  O   ALA A   0      47.947  13.654  24.448  1.00  0.00           O  
ATOM      5  CB  ALA A   0      44.817  13.873  25.699  1.00  0.00           C  
ATOM      6  H1  ALA A   0      45.750  12.989  28.032  1.00  0.00           H  
ATOM      7  H2  ALA A   0      47.376  13.314  27.752  1.00  0.00           H  
ATOM      8  H3  ALA A   0      46.286  12.328  25.836  1.00  0.00           H  
ATOM      9  H4  ALA A   0      47.435  15.250  25.577  1.00  0.00           H  
ATOM     10  H5  ALA A   0      44.230  14.075  26.607  1.00  0.00           H  
ATOM     11  H6  ALA A   0      44.901  14.794  25.104  1.00  0.00           H  
ATOM     12  H7  ALA A   0      44.316  13.096  25.104  1.00  0.00           H  
ATOM     13  N   ALA A   1      50.424  13.298  22.736  1.00  0.00           N  
ATOM     14  CA  ALA A   1      50.417  13.697  21.333  1.00  0.00           C  
ATOM     15  C   ALA A   1      51.805  14.263  21.035  1.00  0.00           C  
ATOM     16  O   ALA A   1      52.408  13.939  20.017  1.00  0.00           O  
ATOM     17  CB  ALA A   1      49.360  14.751  21.095  1.00  0.00           C  
ATOM     18  H1  ALA A   1      50.409  14.138  23.332  1.00  0.00           H  
ATOM     19  H2  ALA A   1      49.593  12.722  22.933  1.00  0.00           H  
ATOM     20  H3  ALA A   1      50.187  12.844  20.678  1.00  0.00           H  
ATOM     21  H4  ALA A   1      52.268  14.968  21.741  1.00  0.00           H  
ATOM     22  H5  ALA A   1      48.988  15.122  22.061  1.00  0.00           H  
ATOM     23  H6  ALA A   1      49.796  15.585  20.525  1.00  0.00           H  
ATOM     24  H7  ALA A   1      48.527  14.313  20.525  1.00  0.00           H  
TER      25      ALA A   1
ENDMDL\n"""
    assert pdb_str == expected
