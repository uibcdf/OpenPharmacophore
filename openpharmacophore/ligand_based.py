from openpharmacophore import Pharmacophore
from openpharmacophore.visualization.view_ligands import view_ligands
from openpharmacophore.extractors.dbscan import dbscan_pharmacophore
from openpharmacophore.io.mol2 import load_mol2_file
from rdkit import Chem
import nglview as nv

class LigandBasedPharmacophore(Pharmacophore):

    def __init__(self, elements=[], molecular_system=None):
        super().__init__(elements=elements, molecular_system=molecular_system)
    
    def from_ligand_list(self, ligands, method, radius=1, feat_list=None, point_type="spheres"):

        """Compute pharmacophore from a list of rdkit molecules 

        Parameters
        ----------
        ligands: :obj: rdkit.Chem.rdmolfiles.SmilesMolSupplier or list of rdkit.Chem.rdchem.Mol
            List of ligands
        
        method: str
            Name of method or algorithm to compute the ligand based pharmacophore

        radius: int (optional)
            Lenght of the radius of the parmacohporic points (Default: 1)
        
        feat_list: list of str (optional)
            List of features that will be used to compute the pharmacophore

        Returns
        -------

        """
        if not isinstance(ligands, list):
            ligands = list(ligands)
        
        if method == "dbscan":
            points, ligands = dbscan_pharmacophore(ligands, radius=radius, feat_list=feat_list)
        else:
            raise NotImplementedError

        self.elements = points
        self.n_elements = len(self.elements)
        self.molecular_system = ligands
        self.extractor = method

    
    def from_ligand_file(self, fname, method, radius=1, feat_list=None, point_type="spheres"):
        """Compute pharmacophore from a file of ligands

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

        Returns
        -------
        """
        accepted_files = ["smi", "mol2", "sdf", "pdb"]
        fextension = fname.split(".")[-1]
        if fextension not in accepted_files:
            raise NotImplementedError
        
        if fextension == "smi":
            ligands = Chem.SmilesMolSupplier(fname, delimiter='\t', titleLine=False)
            ligands = list(ligands)
        elif fextension == "mol2":
            ligands = load_mol2_file(fname)
        elif fextension == "sdf":
            ligands = Chem.SDMolSupplier("../data/abl1/actives_final.sdf")
            ligands = list(ligands)
        elif fextension == "pdb":
            raise NotImplementedError

        self.from_ligand_list(ligands=ligands, method=method, radius=radius, feat_list=feat_list, point_type=point_type)

    def show(self, show_ligands=True, palette="openpharmacophore"):
    
        if self.molecular_system and show_ligands:
            view = view_ligands(self.molecular_system)
        else:
            view = nv.NGLWidget()
        
        self.add_to_NGLView(view, palette=palette)

        return view


