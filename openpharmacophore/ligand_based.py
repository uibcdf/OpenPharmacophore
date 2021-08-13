from openpharmacophore import Pharmacophore
from openpharmacophore.visualization.view_ligands import view_ligands
from openpharmacophore.extractors.dbscan import dbscan_pharmacophore
from openpharmacophore.io.mol2 import load_mol2_file
from rdkit import Chem
import nglview as nv

class LigandBasedPharmacophore(Pharmacophore):

    """ Class to store and compute ligand-based pharmacophores

    Inherits from pharmacophore

    Parameters
    ----------

    elements : :obj:`list` of :obj:`openpharmacophore.pharmacoforic_elements`
        List of pharmacophoric elements

    molecular_system : :obj:`molsysmt.MolSys`
        Molecular system from which this pharmacophore was extracted.

    Attributes
    ----------

    elements : :obj:`list` of :obj:`openpharmacophore.pharmacoforic_elements`
        List of pharmacophoric elements

    n_elements : int
        Number of pharmacophoric elements

    extractor : :obj:`openpharmacophore.extractors`
        Extractor object used to elucidate the pharmacophore

    molecular_system : :obj:`molsysmt.MolSys`
        Molecular system from which this pharmacophore was extracted.

    """

    def __init__(self, elements=[], molecular_system=None):
        super().__init__(elements=elements, molecular_system=molecular_system)
    
    def from_ligand_list(self, ligands, method, radius=1, feat_list=None, feat_def=None, point_type="spheres"):

        """ Compute pharmacophore from a list of rdkit molecules 

        Parameters
        ----------
        ligands: :obj: list of rdkit.Chem.rdchem.Mol rdkit.Chem.SmilesMolSupplier or rdkit.Chem.SDMolSupplier
            List of ligands
        
        method: str
            Name of method or algorithm to compute the ligand based pharmacophore

        radius: float (optional)
            Lenght of the radius of the parmacohporic points (Default: 1)
        
        feat_list: list of str (optional)
            List of features that will be used to compute the pharmacophore
        
        feat_def: dict
            Definitions of the pharmacophoric points. 
            Dictionary which keys are SMARTS strings and values are feature names.

        Note
        -------
        Nothing is returned. The pharmacophore elements are updated with those calculated from the list of ligands.
        The molecular system is updated with the set of ligands and the extractor is updated accoirding to the method used.

        """
        if not isinstance(ligands, list):
            ligands = list(ligands)
        
        if method == "dbscan":
            points, ligands = dbscan_pharmacophore(ligands, radius=radius, feat_list=feat_list, feat_def=feat_def)
        else:
            raise NotImplementedError

        self.elements = points
        self.n_elements = len(self.elements)
        self.molecular_system = ligands
        self.extractor = method

    
    def from_ligand_file(self, fname, method, radius=1, feat_list=None, feat_def=None, point_type="spheres"):

        """ Compute pharmacophore from a file of ligands

        Accepted file formats: smi, mol2, sdf, pdb 

        Parameters
        ----------
        fname: str
            Name of the file containing the ligands
        
        method: str
            Name of method or algorithm to compute the ligand based pharmacophore

        radius: int (optional)
            Lenght of the radius of the parmacohporic points (Default: 1)
        
        feat_list: list of str (optional)
            List of features that will be used to compute the pharmacophore
        
        feat_def: dict
            Definitions of the pharmacophoric points. 
            Dictionary which keys are SMARTS strings and values are feature names.

        Note
        -------
        Nothing is returned. The pharmacophore elements are updated with those calculated from the file of ligands.
        The molecular system is updated with the set of ligands and the extractor is updated accoirding to the method used.

        """
        accepted_files = ["smi", "mol2", "sdf", "pdb"]
        fextension = fname.split(".")[-1]
        if fextension not in accepted_files:
            raise NotImplementedError
        
        if fextension == "smi":
            ligands = Chem.SmilesMolSupplier(fname, delimiter='\t', titleLine=False)
        elif fextension == "mol2":
            ligands = load_mol2_file(fname)
        elif fextension == "sdf":
            ligands = Chem.SDMolSupplier(fname)
        elif fextension == "pdb":
            raise NotImplementedError

        self.from_ligand_list(ligands=ligands, method=method, radius=radius, feat_list=feat_list, feat_def=feat_def, point_type=point_type)

    def show(self, show_ligands=True, palette="openpharmacophore"):
    
        """Showing the pharmacophore model together with the molecular system from with it was
        extracted as a new view (NGLWidget) from NGLView.

        Parameters
        ----------

        show_ligands: bool
            If true the ligands associated to the pharmacophore molecular system are added 
            to the view. (Default: True)

        palette: :obj: `str`, dict
            Color palette name or dictionary. (Default: 'openpharmacophore')

        Returns
        -------
        nglview.NGLWidget
            An nglview.NGLWidget is returned with the 'view' of the pharmacophoric model and the
            molecular system used to elucidate it.
        """

        if self.molecular_system and show_ligands:
            view = view_ligands(self.molecular_system)
        else:
            view = nv.NGLWidget()
        
        self.add_to_NGLView(view, palette=palette)

        return view

