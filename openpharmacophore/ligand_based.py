from openpharmacophore import Pharmacophore
from openpharmacophore.visualization.view_ligands import view_ligands
from openpharmacophore.pharmacophoric_elements.features.color_palettes import get_color_from_palette_for_feature
from openpharmacophore.extractors.dbscan import dbscan_pharmacophore
import nglview as nv
import pyunitwizard as puw 

class LigandBasedPharmacophore(Pharmacophore):

    def __init__(self, elements=[], molecular_system=None):
        super().__init__(elements=elements, molecular_system=molecular_system)
    
    def from_ligand_list(self, ligands, method, radius=1, feat_list=None):
        
        if method == "dbscan":
            tmp_pharmacophore = dbscan_pharmacophore(ligands, radius)
        else:
            raise NotImplementedError

        self.elements = tmp_pharmacophore.elements
        self.n_elements = len(self.elements)
        self.molecular_system = tmp_pharmacophore.molecular_system

    def show(self, show_ligands=True, palette="openpharmacophore"):
    
        if self.molecular_system and show_ligands:
            view = view_ligands(self.molecular_system)
        else:
            view = nv.NGLWidget()
        
        if palette == "openpharmacophore":
            feature_colors = {
                'positive charge': (0.12, 0.36, 0.52), # Blue
                'negative charge': (0.90, 0.30, 0.24),  # Red
                'hb acceptor': (0.90, 0.30, 0.24),  # Red
                'hb donor': (0.13, 0.56, 0.30), # Green
                'included volume': (0, 0, 0), # Black,
                'excluded volume': (0, 0, 0), # Black
                'hydrophobicity': (1, 0.9, 0),  # Yellow
                'aromatic ring': (1, 0.9, 0),  # Yellow
            }
        else:
            raise NotImplementedError

        # TODO: implement visualization of vectors. Openpharmacophore palette is not working

        for i, element in enumerate(self.elements):
            center = puw.get_value(element.center, to_unit="angstroms").tolist()
            radius = puw.get_value(element.radius, to_unit="angstroms")
            # feature_color = get_color_from_palette_for_feature(element.feature_name, color_palette=palette)
            feature_color = feature_colors[element.feature_name]
            label = f"{element.feature_name}_{i}"
            view.shape.add_sphere(center, feature_color, radius, label)
        
        return view


