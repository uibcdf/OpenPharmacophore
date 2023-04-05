# Use this file to define constants used in different modules
from pyunitwizard._private.quantity_or_unit import QuantityLike
from matplotlib.colors import to_rgb


TRAJ_FORMATS = [
    "pdb",
    "h5",
    "xtc",
    "dcd",
    "trr",
    "gro",
    "binpos",
]

TOP_FORMATS = [
    "pdb",
]

MOL_FORMATS = [
    "smi",
    "mol2",
    "sdf",
    "mol",
    "xyz",
]

MOL_STR_FORMATS = [
    "smi",
    "smarts",
    "inchi",
    "mol2",
    "pdb",
]


PALETTE = {
    'positive charge': to_rgb('#3498DB'),  # Blue
    'negative charge': to_rgb('#884EA0'),  # Purple
    'hb acceptor': to_rgb('#B03A2E'),  # Red
    'hb donor': to_rgb('#17A589'),  # Green
    'included volume': to_rgb('#707B7C'),  # Gray
    'excluded volume': to_rgb('#283747'),  # Black
    'hydrophobicity': to_rgb('#F5B041'),  # Orange
    'aromatic ring': to_rgb('#F1C40F'),  # Yellow
}


FEAT_TYPES = frozenset([
    "aromatic ring",
    "hydrophobicity",
    "negative charge",
    "positive charge",
    "hb acceptor",
    "hb donor",
])

FEAT_TO_CHAR = {
    "hb acceptor": "A",
    "hb donor": "D",
    "aromatic ring": "R",
    "hydrophobicity": "H",
    "positive charge": "P",
    "negative charge": "N",
    "excluded volume": "E",
    "included volume": "I",
}

CHAR_TO_FEAT = {
    char: name for name, char in FEAT_TO_CHAR.items()
}
