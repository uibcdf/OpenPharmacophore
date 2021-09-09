from openpharmacophore import Pharmacophore
from openpharmacophore import pharmacophoric_elements as phe
# import molsysmt as msm
import numpy as np
import nglview as nv
from plip.structure.preparation import PDBComplex
import pyunitwizard as puw
from rdkit import Chem

class StructuredBasedPharmacophore(Pharmacophore):

    """ Class to store and compute structured-based pharmacophores

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
    
    @classmethod
    def from_pdb(cls, file_name, radius=1.0, ligand_id=None):
        """Class method to obtain a pharmacophore from a pdb file containing
            a protein-ligand complex.

        Parameters
        ----------
        file_name: str
            Name of the pdb file containing the protein-ligand complex.
        
        radius: float (optional)
            Radius of the spheres of the pharmacophoric points. (Default: 1.0)
        
        ligand_id: str (optional)
            Id of the ligand for which the pharmacophore will be computed

        """
        # TODO: this method should also accept a PDB ID without needing to pass
        # the actual file  
        all_interactions, pdb_string = StructuredBasedPharmacophore._protein_ligand_interactions(file_name)
        if ligand_id:
            for id in all_interactions.keys():
                if ligand_id in id:
                    ligand_id = id
            interactions = all_interactions[ligand_id]
        elif len(all_interactions) == 1:
            interactions = list(all_interactions.values())[0]
        else:
            # If there is more than one ligand and is not specified, prompt the user for the ligand name
            print(f"{file_name} PDB contains the following ligands:\n")
            for lig_id in all_interactions.keys():
                print(lig_id)
            print("\nPlease enter for which one the pharmacophore should be computed ")
            ligand_id = input()
            interactions = all_interactions[ligand_id]
        
        pharmacophoric_points = StructuredBasedPharmacophore._sb_pharmacophore_points(interactions, radius)

        # molecular_system = msm.convert(pdb_string, to_form="molsysmt.MolSys")
        molecular_system = Chem.rdmolfiles.MolFromPDBBlock(pdb_string)

        return cls(elements=pharmacophoric_points, molecular_system=molecular_system)
        
    @staticmethod
    def _protein_ligand_interactions(file_name):
        """Static method to calculate protein-ligand interactions for each ligand in
            the pdb file

        Parameters
        ----------
        file_name: str
            Name of the pdb file containing the protein-ligand complex.

        Returns
        -------
        all_interactions: dict
            Dictionary which keys are ligand Ids and values are all the interaction data 
            for that ligand. 
        
        pdb_string: str
            The corrected pdb stucture as a string.

        """
        mol_system = PDBComplex()
        mol_system.load_pdb(file_name)
        mol_system.analyze()

        pdb_string = mol_system.corrected_pdb
        
        # Dictionary with interactions for each small molecule in the protein
        all_interactions = mol_system.interaction_sets 

        # Filter ions and molecules with less than 5 atoms
        for ligand in mol_system.ligands:
            n_atoms = len(ligand.atomorder)
            if ligand.type != "SMALLMOLECULE" or n_atoms <= 5:
                ligand_id = ":".join([ligand.hetid, ligand.chain, str(ligand.position)])
                all_interactions.pop(ligand_id)
        
        return all_interactions, pdb_string
    
    @staticmethod
    def _sb_pharmacophore_points(interactions, radius):
        """Static method to obtain a list of pharmacophoric points
           for the protein-ligand complex.

        Parameters
        ----------
        interacions: plip.srtucture.preparation.PLInteraction
            object containing all interacion data for a single ligand and a protein.
        
        radius: float (optional)
            Radius of the spheres of the pharmacophoric points. (Default: 1.0)

        Returns
        -------
        points_filtered: list of openpharmacophore.pharmacophoric_elements
            A list of pharmacophoric points.

        """
    
        radius = puw.quantity(radius, "angstroms")

        # List with all interactions
        interactions = interactions.all_itypes
        # list of pharmacophoric_elements
        points = []
        # list oh hydrophobic points
        hydrophobics = [] 
        for interaction in interactions:
            interaction_name = type(interaction).__name__
            
            if interaction_name == "pistack":
                ligand_center = np.array(interaction.ligandring.center)
                protein_center = np.array(interaction.proteinring.center)
                direction = protein_center - ligand_center
                aromatic = phe.AromaticRingSphereAndVector(
                    center=puw.quantity(ligand_center, "angstroms"),
                    radius=radius,
                    direction=direction
                )
                points.append(aromatic)

            elif interaction_name == "hydroph_interaction":
                center = puw.quantity(interaction.ligatom.coords, "angstroms")
                hydrophobic = phe.HydrophobicSphere(
                    center=center,
                    radius=radius
                )
                hydrophobics.append(hydrophobic)
            
            elif interaction_name == "saltbridge":
                if interaction.protispos:
                    # The ligand has a negative charge
                    center = puw.quantity(interaction.negative.center, "angstroms")
                    charge_sphere = phe.NegativeChargeSphere(
                        center=center,
                        radius=radius
                    ) 
                else:
                    # The ligand has a positive charge
                    center = puw.quantity(interaction.positive.center, "angstroms")
                    charge_sphere = phe.PositiveChargeSphere(
                        center=center,
                        radius=radius
                    ) 
                points.append(charge_sphere)
            
            elif interaction_name == "hbond":
                if interaction.protisdon:
                    # The ligand has an acceptor atom
                    ligand_acceptor_center = np.array(interaction.a.coords)
                    protein_donor_center = np.array(interaction.d.coords)
                    direction = ligand_acceptor_center - protein_donor_center 
                    acceptor = phe.HBAcceptorSphereAndVector(
                        center=puw.quantity(ligand_acceptor_center, "angstroms"),
                        radius=radius,
                        direction=direction
                    )
                    points.append(acceptor)
                else:
                    # The ligand has a donor atom
                    ligand_donor_center = np.array(interaction.d.coords)
                    protein_acceptor_center = np.array(interaction.a.coords)
                    direction = protein_acceptor_center - ligand_donor_center
                    donor = phe.HBDonorSphereAndVector(
                        center=puw.quantity(ligand_donor_center, "angstroms"),
                        radius=radius,
                        direction=direction
                    )
                    points.append(donor)
            else:
                # TODO: Incorporate other features such as halogenbonds, waterbridges and metal complexes
                continue

        # Remove duplicate points
        points_filtered = []
        for p in points:
            if p not in points_filtered:
                points_filtered.append(p)
        
        # Group hydrophobic interactions within a distance of 2.5 angstroms
        if hydrophobics:
            skip = []
            for inx, hyd1 in enumerate(hydrophobics):
                if hyd1 in skip:
                    continue
                centroid = hyd1.center
                count = 0
                hyd1_center = puw.get_value(hyd1.center, "angstroms")

                for hyd2 in hydrophobics[inx + 1:]:
                    if hyd1 == hyd2:
                        skip.append(hyd2)
                        continue
                    hyd2_center = puw.get_value(hyd2.center, "angstroms")
                    distance = np.linalg.norm(hyd1_center - hyd2_center)

                    if distance <= 2.5:
                        centroid += hyd2.center
                        count += 1
                        skip.append(hyd2)
                if count > 0:
                    centroid /= (count + 1)
                
                points_filtered.append(phe.HydrophobicSphere(center=centroid, radius=radius))
        
        return points_filtered

    def show(self, palette='openpharmacophore'):

        """Showing the pharmacophore model together with the molecular system from with it was
        extracted as a new view (NGLWidget) from NGLView.

        Parameters
        ----------
        palette: :obj: `str`, dict
            Color palette name or dictionary. (Default: 'openpharmacophore')

        Returns
        -------
        nglview.NGLWidget
            An nglview.NGLWidget is returned with the 'view' of the pharmacophoric model and the
            molecular system used to elucidate it.

        """

        # Molsysmt throws an exception. Temporarily using rdkit 
        # view = msm.view(self.molecular_system, standardize=False)
        view = nv.show_rdkit(self.molecular_system)
        self.add_to_NGLView(view, palette=palette)
        return view
    
