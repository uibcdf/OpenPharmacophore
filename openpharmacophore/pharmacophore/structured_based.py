# OpenPharmacophore
from openpharmacophore._private_tools.exceptions import (FetchError, InvalidFileFormat,
                                                         NoLigandsError, OpenPharmacophoreIOError)
from openpharmacophore import Pharmacophore
from openpharmacophore.io.pharmer import from_pharmer, _pharmer_dict
from openpharmacophore.pharmacophore.color_palettes import get_color_from_palette_for_feature
from openpharmacophore import PharmacophoricPoint
from openpharmacophore.pharmacophore.chemical_features import PharmacophoricPointExtractor
# Third Party
from MDAnalysis.lib.util import NamedStream
import numpy as np
import nglview as nv
from plip.structure.preparation import PDBComplex
import pyunitwizard as puw
from rdkit import Chem, RDLogger
from rdkit.Chem.Draw import rdMolDraw2D
# Standard Library
from collections import defaultdict
import copy
from io import StringIO, BytesIO
import json
import requests
import re
import tempfile
from typing import List, Optional, Tuple
import warnings

RDLogger.DisableLog('rdApp.*')  # Disable rdkit warnings


class StructuredBasedPharmacophore(Pharmacophore):
    """ Class to store and derive pharmacophores from protein-ligand complexes.

    Inherits from pharmacophore

    Parameters
    ----------

    pharmacophoric_points : list of openpharmacophore.PharmacophoricPoint
        List of pharmacophoric pharmacophoric_points

    molecular_system : rdkit.Chem.Mol
        The protein-ligand complex from which this pharmacophore was extracted.

    Attributes
    ----------

    pharmacophoric_points : list of openpharmacophore.PharmacophoricPoint
        List of pharmacophoric pharmacophoric_points

    n_pharmacophoric_points : int
        Number of pharmacophoric pharmacophoric_points

    molecular_system : rdkit.Chem.Mol
        The protein from which this pharmacophore was extracted.
    
    ligand : rdkit.Chem.Mol
        The ligand from the protein-ligand complex.

    """

    def __init__(self, pharmacophoric_points: List[PharmacophoricPoint],
                 molecular_system: Optional[Chem.Mol] = None, ligand: Optional[Chem.Mol] = None) -> None:
        super().__init__(pharmacophoric_points=pharmacophoric_points)
        self.molecular_system = molecular_system
        self.ligand = ligand

    # TODO: add methods to set the molecular system and ligand

    @classmethod
    def from_pdb(cls, pdb: str, radius: float = 1.0, ligand_id: Optional[str] = None,
                 hydrophobics: str = "smarts", load_mol_system: bool = True,
                 load_ligand: bool = True) -> "StructuredBasedPharmacophore":
        """ Class method to obtain a pharmacophore from a pdb file containing a protein-ligand complex. 
            
            Only the interactions of a single ligand will be computed in case  the protein structure contains
            more than one ligand. The id of the ligand can be passed as a string, or if it is not passed a 
            list of  the ligands will appear and the user will be asked to enter which one should be used.  

        Parameters
        ----------
        pdb : str
            PDB id or path to the pdb file containing the protein-ligand complex.
        
        radius : float, optional
            Radius in angstroms of the pharmacophoric points. (Default: 1.0)
        
        ligand_id : str, optional
            Id of the ligand for which the pharmacophore will be computed.
        
        Returns
        -------
        openpharmacophore.StructuredBasedPharmacophore
            The pharmacophore extracted from the protein-ligand complex.

        """
        # TODO: remove hydrophobics argument
        if isinstance(pdb, str):
            pdb_name = pdb
            pattern = re.compile('[0-9][a-zA-Z_0-9]{3}')  # PDB id pattern
            # Check if an ID was passed or a file
            if pdb.endswith(".pdb"):
                as_string = False
            elif pattern.match(pdb):
                pdb = StructuredBasedPharmacophore._fetch_pdb(pdb)
                as_string = True
            else:
                raise OpenPharmacophoreIOError("Invalid file or PDB id")
        # For tempfile.TemporaryFile
        elif isinstance(pdb, tempfile._TemporaryFileWrapper):
            as_string = False
            pdb = pdb.name
        # For mdanalysis streams
        elif isinstance(pdb, NamedStream):
            as_string = True
            pdb = pdb.getvalue()
        else:
            raise TypeError("pdb must be of type str or MDAnalysis.lib.util.NamedStream")

        # pdb_string is the "corrected" pdb that plip generates
        all_interactions, pdb_string, ligands = StructuredBasedPharmacophore._protein_ligand_interactions(pdb,
                                                                                                          as_string=as_string)
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

        if load_ligand:
            ligand_sio = StringIO(ligand_sdf_str)
            ligand_bio = BytesIO(ligand_sio.read().encode("utf8"))
            ligand = [mol for mol in Chem.ForwardSDMolSupplier(ligand_bio)][0]
            if ligand is None:
                warnings.warn("Ligand could not be transformed to rdkit molecule. Using default plip interactions")
                hydrophobics = "plip"
        else:
            ligand = None
            hydrophobics = "plip"

        pharmacophoric_points = StructuredBasedPharmacophore._sb_pharmacophore_points(interactions, radius, ligand,
                                                                                      hydrophobics)

        if load_mol_system:
            molecular_system = Chem.rdmolfiles.MolFromPDBBlock(pdb_string)
        else:
            molecular_system = None

        if ligand:
            ligand = Chem.AddHs(ligand, addCoords=True)
            # If the ligand gets extracted successfully the indices of the pharmacophoric points must be updated to
            # match those of the rdkit molecule. rdkit is zero indexed while pybel indices start at 1. So 1 needs to
            # be subtracted from each index
            for point in pharmacophoric_points:
                indices = {i - 1 for i in point.atom_indices}
                point.atom_indices = indices

        return cls(pharmacophoric_points=pharmacophoric_points,
                   molecular_system=molecular_system, ligand=ligand)

    @classmethod
    def from_file(cls, file_name: str, load_mol_sys: bool = True) -> "StructuredBasedPharmacophore":
        """ Class method to load a structured based pharmacophore from a file.
        
            Currently, supports only json format from pharmer.

        Parameters
        ---------
        file_name : str
            Name of the file containing the pharmacophore.

        load_mol_sys: bool
            Whether to store the molecular system the pharmacophore is associated to

        """
        file_extension = file_name.split(".")[-1]
        if file_extension == "json":
            points, receptor, ligand = from_pharmer(file_name, load_mol_sys)
        else:
            raise InvalidFileFormat(f"Invalid file type, \"{file_name}\" is not a supported file format")

        return cls(points, receptor, ligand)

    @staticmethod
    def _fetch_pdb(pdb_id: str) -> str:
        """ Fetch a protein structure from the Protein Data Bank.
            
            Parameters
            ----------
            pdb_id : str of len 4
                The id of the protein structure.
            
            Returns
            -------
            pdb_str : str
                The pdb as a string.
        """
        url = 'http://files.rcsb.org/download/{}.pdb'.format(pdb_id)
        res = requests.get(url, allow_redirects=True)

        if res.status_code != requests.codes.ok:
            raise FetchError("Could not fetch pdb from {}".format(url))

        pdb_str = res.content.decode()

        return pdb_str

    @staticmethod
    def _protein_ligand_interactions(pdb: str, as_string: bool) -> Tuple[dict, str, dict]:
        """Static method to calculate protein-ligand interactions for each ligand in
            the pdb file.

        Parameters
        ----------
        pdb : str
            File or string containing the protein-ligand complex.

        as_string : bool
            Variable to know if the pdb passed is a string or a file.

        Returns
        -------
        all_interactions : dict
            Dictionary which keys are ligand ids and values are all the interaction data 
            for that ligand. 
        
        pdb_string : str
            The corrected pdb structure as a string.

        ligands : dict
             Dictionary which keys are ligand Ids and values are pybel molecules
        """
        mol_system = PDBComplex()
        mol_system.load_pdb(pdb, as_string=as_string)
        mol_system.analyze()

        pdb_string = mol_system.corrected_pdb

        # Dictionary with interactions for each small molecule in the protein
        all_interactions = mol_system.interaction_sets

        # Store ligands, so we can extract their smiles later on.
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
    def _sb_pharmacophore_points(interactions, radius: float,
                                 ligand: Chem.Mol, hydrophobics: str = "smarts") -> List[PharmacophoricPoint]:
        """Static method to obtain a list of pharmacophoric points for the protein-ligand complex.

        Parameters
        ----------
        interactions : plip.structure.preparation.PLInteraction
            object containing all interaction data for a single ligand and a protein.
        
        radius : float 
            Radius of the spheres of the pharmacophoric points.

        ligand : rdkit.Chem.Mol
            The ligand in from the protein-ligand complex.
        
        hydrophobics : {"plip", "smarts"}
            Can be "plip" or "smarts". If the former is chosen the hydrophobic points will
            be retrieved from the interactions plip calculates. Else smarts patterns will be used. 

        Returns
        -------
        list of openpharmacophore.pharmacophoric_point.PharmacophoricPoint
            A list of pharmacophoric points.

        """
        radius = puw.quantity(radius, "angstroms")

        # List with all interactions
        interactions = interactions.all_itypes
        # list of pharmacophoric_pharmacophoric_points
        points = []
        # list oh hydrophobic points
        hydrophobic_points = []
        for interaction in interactions:
            interaction_name = type(interaction).__name__

            if interaction_name == "pistack":
                ligand_center = np.array(interaction.ligandring.center)
                protein_center = np.array(interaction.proteinring.center)
                direction = protein_center - ligand_center
                atom_indices = {atom.idx for atom in interaction.ligandring.atoms}
                aromatic = PharmacophoricPoint(
                    feat_type="aromatic ring",
                    center=puw.quantity(ligand_center, "angstroms"),
                    radius=radius,
                    direction=direction,
                    atom_indices=atom_indices
                )
                points.append(aromatic)

            elif interaction_name == "hydroph_interaction":
                if hydrophobics != "plip":
                    continue
                center = puw.quantity(interaction.ligatom.coords, "angstroms")
                atom_inx = {interaction.ligatom.idx}
                hydrophobic = PharmacophoricPoint(
                    feat_type="hydrophobicity",
                    center=center,
                    radius=radius,
                    direction=None,
                    atom_indices=atom_inx
                )
                hydrophobic_points.append(hydrophobic)

            elif interaction_name == "saltbridge":
                if interaction.protispos:
                    # The ligand has a negative charge
                    center = puw.quantity(interaction.negative.center, "angstroms")
                    atom_indices = {atom.idx for atom in interaction.negative.atoms}
                    charge_sphere = PharmacophoricPoint(
                        feat_type="negative charge",
                        center=center,
                        radius=radius,
                        atom_indices=atom_indices
                    )
                else:
                    # The ligand has a positive charge
                    center = puw.quantity(interaction.positive.center, "angstroms")
                    atom_indices = {atom.idx for atom in interaction.positive.atoms}
                    charge_sphere = PharmacophoricPoint(
                        feat_type="positive charge",
                        center=center,
                        radius=radius,
                        atom_indices=atom_indices
                    )
                points.append(charge_sphere)

            elif interaction_name == "hbond":
                if interaction.protisdon:
                    # The ligand has an acceptor atom
                    ligand_acceptor_center = np.array(interaction.a.coords)
                    ligand_acceptor_inx = {interaction.a.idx}
                    protein_donor_center = np.array(interaction.d.coords)
                    direction = ligand_acceptor_center - protein_donor_center
                    acceptor = PharmacophoricPoint(
                        feat_type="hb acceptor",
                        center=puw.quantity(ligand_acceptor_center, "angstroms"),
                        radius=radius,
                        direction=direction,
                        atom_indices=ligand_acceptor_inx
                    )
                    points.append(acceptor)
                else:
                    # The ligand has a donor atom
                    ligand_donor_center = np.array(interaction.d.coords)
                    ligand_donor_inx = {interaction.d.idx}
                    protein_acceptor_center = np.array(interaction.a.coords)
                    direction = protein_acceptor_center - ligand_donor_center
                    donor = PharmacophoricPoint(
                        feat_type="hb donor",
                        center=puw.quantity(ligand_donor_center, "angstroms"),
                        radius=radius,
                        direction=direction,
                        atom_indices=ligand_donor_inx
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

        if hydrophobics == "smarts":
            radius = puw.get_value(radius, "angstroms")
            hydrophobic_points = StructuredBasedPharmacophore._smarts_hydrophobics(ligand, radius)

        return points_filtered + hydrophobic_points

    @staticmethod
    def _plip_hydrophobics(hydrophobics: List[PharmacophoricPoint], radius: float) -> List[PharmacophoricPoint]:
        """ Groups plip hydrophobic points that are to close into a single point.

            Parameters
            ----------
            hydrophobics : list of openpharmacophore.PharmacophoricPoint
                List with the hydrophobic points.

            Returns
            -------
            grouped_points : list of openpharmacophore.PharmacophoricPoint
                List with the grouped points.
        """
        # TODO: This method is not working properly
        grouped_points = []
        skip = []
        for inx, hyd1 in enumerate(hydrophobics):
            if hyd1 in skip:
                continue
            centroid = hyd1.center
            indices = hyd1.atom_indices
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
                    indices.union(hyd2.atom_indices)
                    count += 1
                    skip.append(hyd2)
            if count > 0:
                centroid /= (count + 1)

            grouped_points.append(PharmacophoricPoint(
                feat_type="hydrophobicity",
                center=centroid,
                radius=radius,
                direction=None,
                atom_indices=indices))

        return grouped_points

    @staticmethod
    def _smarts_hydrophobics(ligand: Chem.Mol, radius: float) -> List[PharmacophoricPoint]:
        """ Get hydrophobic points usign smarts feature definitions
            
            Parmaeters
            ----------
            ligand : rdkit.Chem.Mol
                The ligand from the protein-ligand complex  

            radius : float
                The radius of the hydrophobic points.

            Returns
            ---------
            points : list of openpharmacophore.PharmacophoricPoint
                List with the hydrophobic points.
        """
        hydrophobic_point_extractor = PharmacophoricPointExtractor(default_radius=radius, features=["hydrophobicity"])
        points = hydrophobic_point_extractor(ligand, 0)

        return points

    def draw(self, file_name: str, img_size: Tuple[int, int] = (500, 500), legend: str = "") -> None:
        """ Draw a 2d representation of the pharmacophore. 
        
            This is a drawing of the ligand with the pharmacophoric features highlighted.

            Parameters
            ----------
            file_name : str
                File where the drawing will be saved. Must be a png file.

            img_size : 2-tuple of int, optional
                The size of the image. (Default=(500,500))

            legend : str, optional
                Image legend.

        """
        if self.ligand is None:
            raise NoLigandsError("Cannot draw pharmacophore if there is no ligand")

        if not file_name.endswith(".png"):
            raise InvalidFileFormat("File must be a png.")

        ligand = copy.deepcopy(self.ligand)
        ligand.RemoveAllConformers()
        ligand = Chem.RemoveHs(ligand)

        atoms = []
        bond_colors = {}
        atom_highlights = defaultdict(list)
        highlight_radius = {}

        for point in self.pharmacophoric_points:

            indices = point.atom_indices
            for idx in indices:

                atoms.append(idx)
                atom_highlights[idx].append(get_color_from_palette_for_feature(point.feature_name))
                highlight_radius[idx] = 0.6

                # Draw aromatic rings bonds
                if point.feature_name == "aromatic ring":
                    for neighbor in ligand.GetAtomWithIdx(idx).GetNeighbors():
                        nbr_idx = neighbor.GetIdx()
                        if nbr_idx not in indices:
                            continue
                        bond = ligand.GetBondBetweenAtoms(idx, nbr_idx).GetIdx()
                        bond_colors[bond] = [get_color_from_palette_for_feature("aromatic ring")]

                # If an atom has more than one feature label will contain both names
                if idx in atoms:
                    if ligand.GetAtomWithIdx(idx).HasProp("atomNote"):
                        label = ligand.GetAtomWithIdx(idx).GetProp("atomNote")
                        label += "|" + str(point.short_name)
                    else:
                        label = point.short_name
                else:
                    label = point.short_name

            ligand.GetAtomWithIdx(idx).SetProp("atomNote", label)

        drawing = rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])
        drawing.DrawMoleculeWithHighlights(ligand, legend, dict(atom_highlights), bond_colors, highlight_radius, {})
        drawing.FinishDrawing()
        drawing.WriteDrawingText(file_name)

    def show(self, palette: str = 'openpharmacophore') -> nv.NGLWidget:
        """ Show the pharmacophore model together with the molecular system from with it was
            extracted.

            Parameters
            ----------
            palette : str, dict
                Color palette name or dictionary. (Default='openpharmacophore')

            Returns
            -------
            nglview.NGLWidget
                An nglview.NGLWidget is returned with the 'view' of the pharmacophoric model and the
                molecular system used to elucidate it.

        """
        view = nv.NGLWidget()
        if self.molecular_system is not None:
            view.add_component(self.molecular_system)
        if self.ligand:
            view.representations = [
                {"type": "cartoon", "params": {
                    "sele": "protein", "color": "residueindex"
                }},
            ]
            view.add_component(self.ligand)
        else:
            view.representations = [
                {"type": "cartoon", "params": {
                    "sele": "protein", "color": "residueindex"
                }},
                {"type": "ball+stick", "params": {
                    "sele": "( not polymer or hetero ) and not ( water or ion )"  # Show ligand
                }}
            ]

        self.add_to_NGLView(view, palette=palette)
        return view

    def to_pharmer(self, file_name: str, save_mol_system: bool = True) -> None:
        """ Save a pharmacophore as a pharmer file (json file). The receptor
            and ligand can also be saved.

            Parameters
            ----------
            file_name: str
                Name of the file that will contain the pharmacophore

            save_mol_system: bool, default=True
                Whether to save the molecular system or just the pharmacophoric points.
    """
        pharmacophore_dict = _pharmer_dict(self.pharmacophoric_points)

        if save_mol_system:
            if self.molecular_system is not None:
                receptor = Chem.MolToPDBBlock(self.molecular_system)
                pharmacophore_dict["receptor"] = receptor

            if self.ligand is not None:
                ligand = Chem.MolToPDBBlock(self.ligand)
                pharmacophore_dict["ligand"] = ligand

        with open(file_name, "w") as outfile:
            json.dump(pharmacophore_dict, outfile)
