from openpharmacophore.screening.screening import VirtualScreening
from rdkit import RDConfig, Chem, Geometry
from rdkit.Chem import ChemicalFeatures, rdDistGeom, rdMolTransforms
from rdkit.Chem.Pharm3D import Pharmacophore, EmbedLib
from rdkit.Numerics import rdAlignment
import pandas as pd
import pyunitwizard as puw
import bisect
import json
from operator import itemgetter
import os
    
class VirtualScreening3D(VirtualScreening):
    """ Class to perform virtual screening with a pharmacophore by 3D alignment of the molecules
        to the pharmacophore. 

    Code adapted from:
    https://github.com/rdkit/UGM_2016/blob/master/Notebooks/Stiefl_RDKitPh4FullPublication.ipynb

    Parameters
    ----------
    pharmacophore: openpharmacophore.Pharmacophore
        The pharmacophore that will be used to screen the database.

    Attributes
    ----------

    aligned_mols : list of 3-tuples (float, rdkit.Chem.mol, str)
        List of molecules that match the pharmacophore. Each tuple is formed by SSD value,
        the molecule object, and the molecule id.

    n_aligned_mols : int
        Number of molecules that match with the pharmacophore. Molecules succesfully aligned to the
        pharmacophore.

    n_fails : int
         Number of molecules that cannot be aligned to the pharmacophore.
    """

    def __init__(self, pharmacophore):
        super().__init__(pharmacophore)
        self.aligned_mols = self.matches 
        self.n_aligned_mols = self.n_matches
        self.scoring_metric = "SSD"
        self._screen_fn = self._align_molecules
        
    def print_report(self):
        """ Prints a summary report of the screening.
        """
        report_str = "Virtual Screening Results\n"
        report_str += "-------------------------\n"
        report_str += "\nMolecules scanned: " 
        report_str += "{:,}".format(self.n_molecules).rjust(36)
        report_str += "\nMolecules matched to pharmacophore: " 
        report_str += str(self.n_aligned_mols).rjust(19)
        report_str += "\nMolecules that didn't match the pharmacophore: " 
        report_str += "{:,}".format(self.n_fails).rjust(8)
        if self.n_aligned_mols > 0:
            report_str += "\nLowest  SSD value: " + str(round(self.aligned_mols[0][0], 4)).rjust(10)
            report_str += "\nHighest SSD value: " + str((round(self.aligned_mols[-1][0], 4))).rjust(10)
            # Calculate mean SSD
            mean = 0
            N = self.n_aligned_mols
            for i in range(N):
                mean += self.aligned_mols[i][0]
            mean /= N
            report_str += "\nAverage SSD value: " + str(round(mean, 4)).rjust(10)         
            # Print top 5 molecules or less if there are less than 5
            if self.n_aligned_mols < 5:
                n_top_mols = min(self.n_aligned_mols, 5)
            else:
                n_top_mols = 5                
            report_str += "\n\nTop {} molecules:\n".format(n_top_mols)
            report_str += "\nZincID " + "SSD".rjust(7)
            report_str += "\n-------  " + " ------\n"
            for i in range(n_top_mols):
                report_str += str(self.aligned_mols[i][-1]) + "   "
                report_str += str(round(self.aligned_mols[i][0], 4)) + "\n"
        print(report_str)

    def save_results_to_file(self, file_name):
        """Save the results of the screening to a file. The file contains the 
           matched molecules ids, smiles and SSD value. File can be saved as 
           csv or json.

           Parameters
           ----------
           file_name: str
                Name of the file

           Notes
           -----
           Does not return anything. A new file is written.
        """
        file_format = file_name.split(".")[-1]

        ssd = [i[0] for i in self.aligned_mols]
        smiles = [Chem.MolToSmiles(i[1]) for i in self.aligned_mols]
        ids = [i[2] for i in self.aligned_mols]

        results = {
            "mol_id": ids,
            "smiles": smiles,
            "ssd": ssd,
        }

        if file_format == "csv":
            df = pd.DataFrame().from_dict(results)
            df.to_csv(file_name, index=False)
        elif file_format == "json":
            json_str = json.dumps(results)
            with open(file_name, "w") as f:
                f.write(json_str)
        else:
            raise NotImplementedError


    def _align_molecules(self, molecules, verbose=0):
        
        """ Align a list of molecules to a given pharmacophore.

        molecules: list of rdkit.Chem.mol
            List of molecules to align.

        verbose: int
            Level of verbosity
        
        Notes
        -------
        Does not return anything. The attributes aligned_mols, n_aligned_mols,
        and n_fails are updated.

        """
        self.n_molecules += len(molecules)

        points = []
        radii = []

        rdkit_element_name = { # dictionary to map openpharmacophore feature names to rdkit feature names
        "aromatic ring": "Aromatic",
        "hydrophobicity": "Hydrophobe",
        "hb acceptor": "Acceptor",
        "hb donor": "Donor",
        "positive charge": "PosIonizable",
        "negative charge": "NegIonizable",
        }

        for element in self.pharmacophore.elements:
            feat_name = rdkit_element_name[element.feature_name]
            center = puw.get_value(element.center, to_unit="angstroms")
            center = Geometry.Point3D(center[0], center[1], center[2])
            points.append(ChemicalFeatures.FreeChemicalFeature(
                feat_name,
                center
            ))
            radius = puw.get_value(element.radius, to_unit="angstroms")
            radii.append(radius)

        rdkit_pharmacophore = Pharmacophore.Pharmacophore(points)
        self._apply_radii_to_bounds(radii, rdkit_pharmacophore)

        fdef = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
        featFactory = ChemicalFeatures.BuildFeatureFactory(fdef)

        for i, mol in enumerate(molecules):

            if verbose == 1 and i % 100 == 0 and i != 0:
                print(f"Screened {i} molecules. Number of matches: {self.n_aligned_mols}; Number of fails: {self.n_fails}")

            bounds_matrix = rdDistGeom.GetMoleculeBoundsMatrix(mol)
            # Check if the molecule features can match with the pharmacophore.
            # TODO: replace this function with a custom one that can take other feature definitions
            can_match, all_matches = EmbedLib.MatchPharmacophoreToMol(mol, featFactory, rdkit_pharmacophore)
            # all_matches is a list of tuples where each tuple contains the chemical features
            if can_match:
                # Match the molecule to the pharmacophore without aligning it
                failed, bounds_matrix_matched, matched_mols, match_details = EmbedLib.MatchPharmacophore(all_matches, 
                                                                                                bounds_matrix,
                                                                                                rdkit_pharmacophore, 
                                                                                                useDownsampling=True)
                if failed:
                    if verbose == 2:
                        print(f"Couldn't embed molecule {i}")
                    self.n_fails += 1
                    continue
            else:
                if verbose == 2:
                    print(f"Couldn't match molecule {i}")
                self.n_fails += 1
                continue
            atom_match = [list(x.GetAtomIds()) for x in matched_mols]
            try:
                mol_H = Chem.AddHs(mol)
                # Embed molecule onto the pharmacophore
                # embeddings is a list of molecules with a single conformer
                b_matrix, embeddings, num_fail = EmbedLib.EmbedPharmacophore(mol_H, atom_match, rdkit_pharmacophore, count=10)
            except Exception as e:
                if verbose == 2:
                    print(e)
                    print (f"Bounds smoothing failed for molecule {i}")
                self.n_fails += 1
                continue
            # Align embeddings to the pharmacophore 
            SSDs = self._transform_embeddings(rdkit_pharmacophore, embeddings, atom_match) 
            best_fit_index = min(enumerate(SSDs), key=itemgetter(1))[0]

            mol_id = mol.GetProp("_Name")
            matched_mol = (SSDs[best_fit_index], embeddings[best_fit_index], mol_id)
            # Append to list in ordered manner
            bisect.insort(self.aligned_mols, matched_mol) 
            self.n_aligned_mols += 1
    
    def _apply_radii_to_bounds(self, radii, pharmacophore):
        """
            Apply the radius of the pharmacophoric points to the bound
            matrix of the pharmacophore

            Parameters
            ----------
            radii: list of float
                list with the radius of each of the pharmacophoric points
            
            pharmacophore: rdkit.Chem.Pharm3D.Pharmacophore
                The pharmacophore wich bounds matriz will be modified
            
            Returns
            -------
            Nothing is returned. The pharmacophore bounds matrix is updated

        """
        for i in range(len(radii)):
            for j in range(i + 1, len(radii)):
                sum_radii = radii[i] + radii[j]
                pharmacophore.setLowerBound(i, j, max(pharmacophore.getLowerBound(i, j) - sum_radii, 0))
                pharmacophore.setUpperBound(i, j, pharmacophore.getUpperBound(i, j) + sum_radii)

    def _get_transform_matrix(self, align_ref, conf_embed, atom_match):

        """Get the transformation matrix for a conformer.

            Parameters
            ----------
            align_ref: list of rdkit.Geometry.Point3D
                list of pharmacophore reference points for the alignment.

            conf_embed: rdkit.Chem.Conformer
                The conformer embedding.

            atom_match: list of list
                List of list of atoms ids that match the pharmacophore.

            Returns
            -------
            a 2-tuple

            ssd: float
                SSD value for the alignment.

            transform_matrix: numpy.ndarray; shape(4, 4)
                The transform matrix.

            """
    
        align_probe = [] # probe points to align to reference points
        for match_ids in atom_match:
            dummy_point = Geometry.Point3D(0.0,0.0,0.0)
            for id in match_ids:
                dummy_point += conf_embed.GetAtomPosition(id)
            dummy_point /= len(match_ids)
            align_probe.append(dummy_point)
        
        ssd, transform_matrix = rdAlignment.GetAlignmentTransform(align_ref, align_probe)
        return ssd, transform_matrix

    def _transform_embeddings(self, pharmacophore, embeddings, atom_match):

        """Transform embeddings. Performs the alignment of the molecules 
            to the pharmacophore.

            Parameters
            ----------
            pharmacophore: rdkit.Chem.Pharm3D.Pharmacophore
                A pharmacophore object.

            embeddings: list of rdkit.Chem.Mol
                List of molecules with a single conformer.

            atom_match: list of list
                List of list of atoms ids that match the pharmacophore.

            Returns
            -------
            SSDs: list of float
                List of sum of square deviations (SSD) values for the alignments.

            """

        align_ref = [f.GetPos() for f in pharmacophore.getFeatures()]
        ssds = []
        for embedding in embeddings:
            conformer = embedding.GetConformer()
            ssd, transform_matrix = self._get_transform_matrix(align_ref, conformer, atom_match)
            # Transform the coordinates of the conformer
            rdMolTransforms.TransformConformer(conformer, transform_matrix)
            ssds.append(ssd)
        
        return ssds


    
