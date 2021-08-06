from rdkit import Chem
import nglview as nv

def view_conformers(molecule):
    """
    Generate a view of the conformers of a molecule.

    Parameters
    -----------
    molecule: an rdkit.Chem.rdchem.Mol

    Returns
    ----------
    view: an nglview.widget.NGLWidget
    """
    view = nv.NGLWidget()
    for conformer in range(molecule.GetNumConformers()):
        mol_string = Chem.MolToMolBlock(molecule, confId=conformer)
        temp_mol = Chem.MolFromMolBlock(mol_string, removeHs=False)
        component = view.add_component(temp_mol)
        component.clear()
        component.add_ball_and_stick(multipleBond=True)
    return view