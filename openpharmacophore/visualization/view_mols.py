import nglview as nv
from rdkit import Chem
from typing import TypeVar, List

Molecule = TypeVar("Molecule", List[Chem.Mol], Chem.Mol)

def view_ligands(molecules: Molecule) -> nv.NGLWidget:
    """
    Generate a view of a set of ligands

    Parameters
    -----------
    molecules : rdkit.Chem.Mol or list of rdkit.Chem.rdchem.Mol
        The molecules that will be visualized.

    Returns
    ----------
    view: nglview.widget.NGLWidget
    """
    if not isinstance(molecules, list):
        molecules = [molecules]
    
    view = nv.NGLWidget()
    for molecule in molecules:
        component = view.add_component(molecule)
        component.clear()
        component.add_ball_and_stick(multipleBond=True)
    return view


def view_conformers(molecule: Chem.Mol) -> nv.NGLWidget:
    """
    Generate a view of the conformers of a molecule.

    Parameters
    -----------
    molecule : rdkit.Chem.Mol
        The molecule which conformers will be visualized.

    Returns
    ----------
    view : nglview.widget.NGLWidget
    """
    view = nv.NGLWidget()
    for conformer in range(molecule.GetNumConformers()):
        mol_string = Chem.MolToMolBlock(molecule, confId=conformer)
        temp_mol = Chem.MolFromMolBlock(mol_string, removeHs=False)
        component = view.add_component(temp_mol)
        component.clear()
        component.add_ball_and_stick(multipleBond=True)
    return view