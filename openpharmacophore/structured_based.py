from openpharmacophore._private_tools.exceptions import FetchError
from openpharmacophore import Pharmacophore
from openpharmacophore import pharmacophoric_elements as phe
from openpharmacophore.utils import ligand_features
# import molsysmt as msm
import numpy as np
import nglview as nv
from plip.structure.preparation import PDBComplex
import pyunitwizard as puw
from rdkit import Chem, RDLogger
from io import StringIO, BytesIO
import requests
import re
import warnings

RDLogger.DisableLog('rdApp.*') # Disable rdkit warnings

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

    def __init__(self, elements=[], molecular_system=None, ligand=None):
        super().__init__(elements=elements, molecular_system=molecular_system)
        self.ligand = ligand
    
    @classmethod
    def from_pdb(cls, pdb, radius=1.0, ligand_id=None, hydrophobics="rdkit"):
        """ Class method to obtain a pharmacophore from a pdb file containing
            a protein-ligand complex. 
            
            Only the interactions of a single lignad will be computed in case 
            the protein structure contains more than one ligand. The id of the 
            ligand can be passed as a string, or if it is not passed a list of 
            the ligands will appear and the user will be asked to enter whicho
            one shoul be used.  

        Parameters
        ----------
        pdb: str
            PDB id or path to the pdb file containing the protein-ligand complex.
        
        radius: float
            Radius of the spheres of the pharmacophoric points. (Default: 1.0)
        
        ligand_id: str (optional)
            Id of the ligand for which the pharmacophore will be computed.
        
        Returns
        -------
        An openpharmacophore.StructuredBasedPharmacophore with the elements 

        """
        pdb_name = pdb
        pattern= re.compile('[0-9][a-zA-Z_0-9]{3}') #PDB id pattern
        # Check if an ID was passed or a file
        if pdb.endswith(".pdb"):
            as_string = False
        elif pattern.match(pdb):
            pdb = StructuredBasedPharmacophore._fetch_pdb(pdb)
            as_string = True
        else:
            raise Exception("Invalid file or PDB id")

        # pdb_string is the "corrected" pdb that plip generates
        all_interactions, pdb_string, ligands = StructuredBasedPharmacophore._protein_ligand_interactions(pdb, as_string=as_string)
        if ligand_id:
            interactions = all_interactions[ligand_id]
            # Extract the ligand and convert to sdf to get its 3D coordinates
            ligand_sdf_str = ligands[ligand_id].write("sdf")
            
        
        elif len(all_interactions) == 1:
            interactions = list(all_interactions.values())[0]
            ligand_pybel = list(ligands.values())[0]
            ligand_sdf_str = ligand_pybel.write("sdf")
        
        else:
            # If there is more than one ligand and is not specified, prompt the user for the ligand name
            print(f"{pdb_name} PDB contains the following ligands:\n")
            for lig_id in all_interactions.keys():
                print(lig_id)
            print("\nPlease enter for which one the pharmacophore should be computed ")
            ligand_id = input()
            interactions = all_interactions[ligand_id]
            ligand_sdf_str = ligands[ligand_id].write("sdf")
        
        ligand_sio = StringIO(ligand_sdf_str)
        ligand_bio = BytesIO(ligand_sio.read().encode("utf8"))
        ligand = [mol for mol in Chem.ForwardSDMolSupplier(ligand_bio)][0]
        if ligand is None:
            warnings.warn("Ligand could not be transformed to rdkit molecule. Using default plip interactions")
            hydrophobics="plip"
        pharmacophoric_points = StructuredBasedPharmacophore._sb_pharmacophore_points(interactions, radius, ligand, hydrophobics)

        # molecular_system = msm.convert(pdb_string, to_form="molsysmt.MolSys")
        molecular_system = Chem.rdmolfiles.MolFromPDBBlock(pdb_string)

        return cls(elements=pharmacophoric_points, molecular_system=molecular_system, ligand=ligand)
    
    @staticmethod
    def _fetch_pdb(pdb_id):
        """ Fetch a protein estructure from PDB.
            
            Parameters
            ----------
            pdb_id: str of len 4
                The id of the protein structure.
            
            Returns
            -------
            pdb_str: str
                The pdb as a string.
        """
        url = 'http://files.rcsb.org/download/{}.pdb'.format(pdb_id)
        res = requests.get(url, allow_redirects=True)
     
        if res.status_code != requests.codes.ok:
            raise FetchError("Could not fetch pdb from {}".format(url)) 
        
        pdb_str = res.content.decode()

        return pdb_str
        
    @staticmethod
    def _protein_ligand_interactions(pdb, as_string):
        """Static method to calculate protein-ligand interactions for each ligand in
            the pdb file

        Parameters
        ----------
        pdb: str
            File or string containing the protein-ligand complex.

        as_string: bool
            Variable to know if the pdb passed is a string or a file.

        Returns
        -------
        all_interactions: dict
            Dictionary which keys are ligand Ids and values are all the interaction data 
            for that ligand. 
        
        pdb_string: str
            The corrected pdb stucture as a string.

        ligands: dict
             Dictionary which keys are ligand Ids and values are pybel molecules
        """
        mol_system = PDBComplex()
        mol_system.load_pdb(pdb, as_string=as_string)
        mol_system.analyze()

        pdb_string = mol_system.corrected_pdb
        
        # Dictionary with interactions for each small molecule in the protein
        all_interactions = mol_system.interaction_sets 

        # Store ligands so we can extract their smiles later on.
        ligands = dict()
      
        for ligand in mol_system.ligands:
            n_atoms = len(ligand.atomorder)
            ligand_id = ":".join([ligand.hetid, ligand.chain, str(ligand.position)])
            # Filter ions and molecules with less than 5 atoms
            if ligand.type != "SMALLMOLECULE" or n_atoms <= 5:
                all_interactions.pop(ligand_id)
                continue
            # ligand.mol is an openbabel molecule
            ligands[ligand_id] = ligand.mol
        
        return all_interactions, pdb_string, ligands
    
    @staticmethod
    def _sb_pharmacophore_points(interactions, radius, ligand, hydrophobics="rdkit"):
        """Static method to obtain a list of pharmacophoric points
           for the protein-ligand complex.

        Parameters
        ----------
        interacions: plip.srtucture.preparation.PLInteraction
            object containing all interacion data for a single ligand and a protein.
        
        radius: float 
            Radius of the spheres of the pharmacophoric points.

        ligand: rdkit.Chem.mol
            The ligand in from the protein-ligand complex.
        
        hydrophobics: str
            Can be "plip" or "rdkit". If the former is chosen the hydrophobic points will
            be retrieved from the interactions plip calculates. Else smarts patterns from
            rdkit will be used. 

        Returns
        -------
        list of openpharmacophore.pharmacophoric_elements
            A list of pharmacophoric points.

        """   
        radius = puw.quantity(radius, "angstroms")

        # List with all interactions
        interactions = interactions.all_itypes
        # list of pharmacophoric_elements
        points = []
        # list oh hydrophobic points
        hydrophobic_points = [] 
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
                if hydrophobics != "plip":
                    continue
                center = puw.quantity(interaction.ligatom.coords, "angstroms")
                hydrophobic = phe.HydrophobicSphere(
                    center=center,
                    radius=radius
                )
                hydrophobic_points.append(hydrophobic)
            
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
        if len(hydrophobic_points) > 1:
            hydrophobic_points = StructuredBasedPharmacophore._plip_hydrophobics(hydrophobic_points, radius)
       
        if hydrophobics == "rdkit":
            radius = puw.get_value(radius, "angstroms")
            hydrophobic_points = StructuredBasedPharmacophore._rdkit_hydrophobics(ligand, radius)
        
        return points_filtered + hydrophobic_points

    @staticmethod
    def _plip_hydrophobics(hydrophobics, radius):
        """ Groups plip hydrophobic points that are to close into a single
            point.

            Parameters
            ----------
            hydrophobics: list of openpharmacophore.pharmacophoric_elements.HydrophobicSphere
                List with the hydrophobic points.

            Returns
            -------
            grouped_points: list of openpharmacophore.pharmacophoric_elements.HydrophobicSphere
                List with the grouped points.
        """
        grouped_points = []
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
            
            grouped_points.append(phe.HydrophobicSphere(center=centroid, radius=radius))
        
        return grouped_points

    @staticmethod
    def _rdkit_hydrophobics(ligand, radius):
        """ Get hydrophobic points usign rdkit feature definitions
            
            Parmaeters
            ----------
            ligand: rdkit.Chem.mol
                The ligand from the protein-ligand complex  

            radius: float
                The radius of the hydrophobic points.

            Returns
            ---------
            points: list of openpharmacophore.pharmacophoric_elements.HydrophobicSphere
                List with the hydrophobic points.
        """
        hydrophobic_smarts = [
        "[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]",
        "[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
        "[C&r3]1~[C&r3]~[C&r3]1",
        "[C&r4]1~[C&r4]~[C&r4]~[C&r4]1",
        "[C&r5]1~[C&r5]~[C&r5]~[C&r5]~[C&r5]1",
        "[C&r6]1~[C&r6]~[C&r6]~[C&r6]~[C&r6]~[C&r6]1",
        "[C&r7]1~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]1",
        "[C&r8]1~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]1",
        "[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
        "[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
        "[$([CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
        "[$([S]~[#6])&!$(S~[!#6])]",
        ]

        smarts_dict = {smarts : "Hydrophobe" for smarts in hydrophobic_smarts}
        points_dict = ligand_features.ligands_pharmacophoric_points(ligand, radius, feat_list=None, feat_def=smarts_dict)
        points = points_dict["ligand_0"]["conformer_0"]
        return points 

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
    
