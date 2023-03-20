from copy import deepcopy
import math
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG

from openpharmacophore.constants import PALETTE


def draw_ligands(ligands, n_per_row, pretty=True):
    """ Draw multiple ligands

        Parameters
        ----------
        ligands : list[Ligand]
            A list with ligands.
        n_per_row : int
            Number of ligands to draw on each row.

        pretty : bool, default=True
            If true the ligands is drawn with no conformers.
    """
    if pretty:
        # Copy to avoid removing conformers from the original molecules
        mols = [deepcopy(lig.to_rdkit()) for lig in ligands]
        for mol in mols:
            mol.RemoveAllConformers()
    else:
        mols = [lig.to_rdkit() for lig in ligands]
    return MolsToGridImage(mols, molsPerRow=n_per_row)


def _atom_highlights(feat_ind):
    """ Returns the atoms that will be highlighted in a drawing
        along with their color and the radius of the circle.

        Parameters
        ----------
        feat_ind : list[dict[str, list[tuple[int]]]]
            Indices of each chemical feature

        Returns
        -------
        atoms : list[tuple[int]]
        colors : list[dict[int, tuple[int]]
        radii : list[dict[int, float]]
    """
    atoms = []
    colors = []
    radii = []

    for lig_feats in feat_ind:
        lig_ind = []
        lig_color = {}
        lig_radii = {}
        for feat_type, ind_list in lig_feats.items():
            color = PALETTE[feat_type]
            for ind_tup in ind_list:
                for ind in ind_tup:
                    lig_ind.append(ind)
                    lig_color[ind] = color
                    lig_radii[ind] = 0.5
        atoms.append(tuple(lig_ind))
        colors.append(lig_color)
        radii.append(lig_radii)

    return atoms, colors, radii


def _drawing_size(mol_width, mol_height, n_mols):
    """ Returns the width and height of a drawing consisting of n_mols
        with the given width and height

       Parameters
       ----------
       mol_width : int
       mol_height : int
       n_mols : int

       Returns
       -------
       width : int
       height : int
           """
    max_per_row = 4
    max_width = max_per_row * mol_width
    if mol_width * n_mols <= max_width:
        width = mol_width * n_mols
        n_rows = 1
    else:
        width = max_width
        n_rows = math.ceil(mol_width * n_mols / width)
    height = n_rows * mol_height
    return width, height


def draw_ligands_chem_feats(ligands, lig_size):
    """ Draw the ligands with their chemical features highlighted.

        Parameters
        ----------
        ligands : list[Ligand]
            A list of ligands.

        lig_size : tuple[int, int]
            A tuple with the width and height of each ligand drawing.

        Returns
        -------
        drawing
    """
    molecules = [deepcopy(lig.to_rdkit()) for lig in ligands]
    for mol in molecules:
        mol.RemoveAllConformers()

    feat_ind = []
    for lig in ligands:
        feat_ind.append(lig.feat_ind)

    atoms, colors, radii = _atom_highlights(feat_ind)
    width, height = _drawing_size(lig_size[0], lig_size[1], len(ligands))
    drawing = MolDraw2DSVG(width=width, height=height,
                           panelWidth=lig_size[0], panelHeight=lig_size[1])
    drawing.DrawMolecules(molecules,
                          highlightAtoms=atoms,
                          highlightAtomColors=colors,
                          highlightAtomRadii=radii
                          )
    drawing.FinishDrawing()
    return drawing
