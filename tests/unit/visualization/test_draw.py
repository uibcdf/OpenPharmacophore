from openpharmacophore.constants import PALETTE
from openpharmacophore.visualization.draw import _atom_highlights, _drawing_size


def test_atom_highlights():
    feat_ind = [
        {"aromatic ring": [(1, 2, 3)]},
        {"hb donor": [(1,)]},
    ]
    ring_color = PALETTE["aromatic ring"]
    donor_color = PALETTE["hb donor"]
    atoms, colors, radii = _atom_highlights(feat_ind)

    expected_atoms = [(1, 2, 3), (1,)]
    expected_colors = [
        {1: ring_color, 2: ring_color, 3: ring_color},
        {1: donor_color}
    ]
    expected_radii = [
        {1: 0.5, 2: 0.5, 3: 0.5},
        {1: 0.5}
    ]
    assert atoms == expected_atoms
    assert expected_colors == colors
    assert radii == expected_radii


def test_drawing_size():
    assert _drawing_size(300, 280, 2) == (600, 280)
    assert _drawing_size(300, 280, 4) == (1200, 280)
    assert _drawing_size(300, 280, 6) == (1200, 560)
