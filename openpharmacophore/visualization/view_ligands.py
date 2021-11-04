import nglview as nv

def view_ligands(molecules):
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