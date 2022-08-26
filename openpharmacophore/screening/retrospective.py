# OpenPharmacophore
from openpharmacophore import Pharmacophore, LigandBasedPharmacophore, StructuredBasedPharmacophore
from openpharmacophore.screening.alignment import apply_radii_to_bounds, transform_embeddings
from openpharmacophore._private_tools.exceptions import BadShapeError, MissingParameters, OpenPharmacophoreValueError
from openpharmacophore._private_tools.screening_arguments import check_virtual_screening_kwargs, is_3d_pharmacophore
# Third Party
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem, RDConfig, DataStructs
from rdkit.Chem import ChemicalFeatures, rdDistGeom
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D
from rdkit.Chem.Pharm2D.Generate import Gen2DFingerprint
from rdkit.Chem.Pharm3D import EmbedLib
from tqdm.auto import tqdm
# Standard library
from collections import namedtuple
from operator import itemgetter
import os
from typing import Tuple, List, Optional, TypeVar

PharmacophoreType = TypeVar("PharmacophoreType", LigandBasedPharmacophore,
                            StructuredBasedPharmacophore, Pharmacophore, DataStructs.SparseBitVect)

MolScore = namedtuple("MolScore", ["score", "id", "mol"])


class RetrospectiveScreening:
    """ Class for performing retrospective virtual screening. 
    
        This class expects molecules classified as actives and inactives. 
        With this class pharmacophore models can be validated.

    Parameters
    ----------
    pharmacophore : Pharmacophore
        The pharmacophore that will be used to screen the database. Can be a Pharmacophore, 
        StructuredBasedPharmacophore, LigandBasedPharmacophore or a fingerprint

    Attributes
    ----------
    matches : list of 3-tuples (float, str, rdkit.Chem.Mol)
        List of molecules that match the pharmacophore. Each tuple is formed by scoring 
        value, the molecule id, and the molecule object.
    
    n_actives : int
        Number of active molecules.

    n_inactives: int
        Number of inactives.

    n_molecules: int
        Number of molecules screened.
    
    n_fails : int
        Number of molecules that cannot be matched to the pharmacophore.
    
    scoring_metric : str
        Metric used to score the molecules, how well they fit to the pharmacophore.
    
    pharmacophore : Pharmacophore
        The pharmacophore that will be used to screen the database. Can be a Pharmacophore, 
        StructuredBasedPharmacophore, LigandBasedPharmacophore or a fingerprint



    """

    def __init__(self, pharmacophore: PharmacophoreType, **kwargs) -> None:

        if is_3d_pharmacophore(pharmacophore):
            self.scoring_metric = "SSD"
            self._screen_fn = self._align_molecules
        elif isinstance(pharmacophore, DataStructs.SparseBitVect):  # For pharmacophore fingerprints
            self.scoring_metric = "Similarity"
            self.similarity_fn, _ = check_virtual_screening_kwargs(**kwargs)
            self._factory = Gobbi_Pharm2D.factory
            self._screen_fn = self._fingerprint_similarity
        else:
            raise TypeError("pharmacophore must be of type Pharmacophore, StructuredBasedPharmacophore, "
                            "LigandBasedPharmacophore, or rdkit.DataStructs.SparseBitVect")

        self.n_actives = 0
        self.n_inactives = 0
        self.bioactivities = None
        self.molecules = []
        self.n_molecules = 0
        self.pharmacophore = pharmacophore

    def from_bioactivity_data(self, smiles: List[Tuple[int, str]], activity: np.ndarray) -> None:
        """ Retrospective screening from a set of molecules classified as active or inactive.
        
            Parameters
            ----------
            smiles : List of 2-tuples
                A list with the molecules for screening. Each element of the list is 
                a tuple, where the first elements is the compound id and the second 
                the smiles of the molecule.

            activity : numpy.ndarray 
                Array with the labels of each molecule; 1 corresponds to an active molecule
                and 0 to an inactive one. An array of rank 1 where the first dimension is 
                equal to the length of the smiles list is expected. 
        
        """
        if len(activity.shape) > 1:
            raise BadShapeError("activity must be an array of rank 1")
        if len(smiles) != activity.shape[0]:
            raise OpenPharmacophoreValueError("smiles and activity must contain the same number of entries")

        self.n_actives = np.sum(activity)
        self.n_inactives = activity.shape[0] - self.n_actives
        self.bioactivities = activity

        molecules = []
        for id_, smi in smiles:
            mol = Chem.MolFromSmiles(smi)
            mol.SetProp("_Name", str(id_))
            molecules.append(mol)

        if self.scoring_metric == "Similarity":
            self._fingerprint_similarity(molecules)
        elif self.scoring_metric == "SSD":
            self._align_molecules(molecules)
        else:
            raise NotImplementedError

    def confusion_matrix(self, threshold: Optional[float] = None) -> np.ndarray:
        """ Compute a confusion matrix
        
            Parameters
            ----------
            threshold : float, optional
                The scoring value from which a molecule will be considered as active.
                Required if the screening was done with fingerprints.
            
            Returns
            -------
            cf_matrix: np.ndarray of shape (2, 2)
                The confusion matrix.

        """
        if self.scoring_metric == "Similarity" and threshold is None:
            raise MissingParameters("Expected a threshold value.")

        if self.scoring_metric == "SSD":
            # TODO: is threshold value of 0 correct? An SSD of 0 means that
            #  the pharmacophore aligns perfectly with the ligand. However this is really
            #  unlikely. We should establish a higher threshold
            threshold = 0.0

        true_positives = 0
        true_negatives = 0
        false_positives = 0
        false_negatives = 0

        for ii, mol in enumerate(self.molecules):
            if self.bioactivities[ii] == 1 and mol.score > threshold:
                true_positives += 1
            elif self.bioactivities[ii] == 1 and mol.score <= threshold:
                false_negatives += 1
            elif self.bioactivities[ii] == 0 and mol.score > threshold:
                false_positives += 1
            elif self.bioactivities[ii] == 0 and mol.score <= threshold:
                true_negatives += 1

        cf_matrix = np.array([[true_positives, false_positives],
                              [false_negatives, true_negatives]])

        assert np.sum(cf_matrix, axis=None) == self.n_molecules

        return cf_matrix

    def auc(self) -> float:
        """ Calculate ROC area under the curve.
        
            Returns
            -------
            float
                The value of the area under the curve.

            References
            ----------
            Fawcett, T. An Introduction to ROC Analysis. Pattern Recognition Letters 2006, 27, 861−874
        """
        scores = [x[0] for x in self.molecules]
        scores = np.array(scores)

        return self._get_auc(scores, self.bioactivities)

    def roc_plot(self, ax: Optional[plt.Axes] = None, label: str = "",
                 random_line: bool = True) -> plt.Axes:
        """ Plot the ROC curve. 
        
            Parameters
            ----------
            ax : matplotlib.axes._subplots.AxesSubplot, optional (Default = None)
                An axes object where the plot will be drawn.

            random_line : bool, default=True
                Whether to plot the line corresponding to a random classifier.

            label : str, default=""
                The label of the ROC curve

            Returns
            ----------
            ax : matplotlib.axes._subplots.AxesSubplot, optional (Default = None)
                An axes object with the plot.

            References
            ----------
            Fawcett, T. An Introduction to ROC Analysis. Pattern Recognition Letters 2006, 27, 861−874

        """
        scores = [x[0] for x in self.molecules]
        scores = np.array(scores)
        fpr, tpr = self._roc_points(scores, self.bioactivities)

        # Plot the curve

        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(fpr, tpr, label=label)
        if random_line:
            ax.plot([0, 1], [0, 1], color="black", linestyle="dashed", label="Random")

        ax.set_ylabel("Sensitivity")
        ax.set_xlabel("1 - Specificity")
        if label or random_line:
            ax.legend()

        return ax

    def enrichment_factor(self, percentage: float) -> float:
        """ Get enrichment factor for the x% of the screened database 

            Parameters
            ----------
            percentage : float
                Percentage of the screened database. Must be between 0 and 100
            
            Returns
            -------
            float
                The enrichment factor
        """
        if percentage < 0 or percentage > 100:
            raise OpenPharmacophoreValueError("percentage must be a number between 0 and 100")

        scores = [x[0] for x in self.molecules]
        scores = np.array(scores)

        return self._calculate_enrichment_factor(scores, self.bioactivities, percentage)

    def ideal_enrichment_factor(self, percentage: float) -> float:
        """ Calculate ideal enrichment factor for the x% of the screened database 

            Parameters
            ----------
            percentage : float
                Percentage of the screened database. Must be between 0 and 100
            
            Returns
            -------
            float
                The ideal enrichment factor"""

        percentage = percentage / 100

        n_molecules = self.bioactivities.shape[0]
        ratio_actives = self.n_actives / n_molecules
        if percentage <= ratio_actives:
            return (100 / ratio_actives) * percentage
        else:
            return 100.0

    def enrichment_plot(self, ax: Optional[plt.Axes] = None, label: str = "",
                        random_line: bool = True, ideal: bool = False) -> plt.Axes:
        """ Create an enrichment plot 
            
            Parameters
            ----------
            ax : matplotlib.axes._subplots.AxesSubplot, optional (Default = None)
                An axes object where the plot will be drawn.

            random_line : bool, default=True
                Whether to plot the line corresponding to a random classifier.

            ideal : bool, default=False
                Whether to plot the ideal enrichment curve

             label : str, default=""
                The label of the enrichment curve
            
            Returns
            ----------
            ax : matplotlib.axes._subplots.AxesSubplot, optional (Default = None)
                An axes object with the plot.
        
        """
        scores = [x[0] for x in self.molecules]
        scores = np.array(scores)

        screened_percentage, percentage_actives_found = self._enrichment_data(scores, self.bioactivities)

        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(screened_percentage, percentage_actives_found, label=label)
        if random_line:
            ax.plot([0, 1], [0, 1], color="black", linestyle="dashed", label="Random")
        if ideal:
            n_molecules = self.bioactivities.shape[0]
            ratio_actives = self.n_actives / n_molecules
            ax.plot([0, ratio_actives, 1], [0, 1, 1], color="red", linestyle="dashed", label="Ideal")

        ax.set_xlabel("% Database Screened")
        ax.set_ylabel("% Actives")
        ax.legend()

        return ax

    @staticmethod
    def _sort_scores_and_labels_for_roc(scores: np.ndarray,
                                        labels: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """ Sorts scores and labels needed for ROC an AUC calculations.
        """
        if len(scores.shape) > 1 or len(labels.shape) > 1:
            raise BadShapeError("scores and labels must be an array of rank one.")
        if scores.shape[0] != labels.shape[0]:
            raise BadShapeError("scores and labels must have the same shape.")
        # Sort scores in descending order and then sort labels
        indices = np.argsort(scores)[::-1]
        scores = np.sort(scores)[::-1]
        labels = labels[indices]

        return scores, labels

    @staticmethod
    def _roc_points(scores: np.ndarray, labels: np.ndarray) -> Tuple[List[float], List[float]]:
        """ Calculate points to plot an ROC curve.
        
            Parameters
            ----------
            scores : np.ndarray of shape(n_samples,)
                The score of each sample.
            
            labels : np.ndarray of shape(n_samples,)
                The ground truth.
            
            Returns
            -------
            false_positive_rate : list of float
                A list with the values of false positive rate.
            
            true_positive_rate : list of float
                A list with the values of true positive rate.
        """
        scores, labels = RetrospectiveScreening._sort_scores_and_labels_for_roc(scores, labels)

        n_positives = np.sum(labels)
        n_negatives = labels.shape[0] - n_positives

        score_prev = -10000000

        false_positives = 0
        true_positives = 0
        false_positive_rate = []
        true_positive_rate = []

        # Calculate points for the ROC plot
        ii = 0
        while ii < labels.shape[0]:

            if scores[ii] != score_prev:
                false_positive_rate.append(false_positives / n_negatives)
                true_positive_rate.append(true_positives / n_positives)
                score_prev = scores[ii]

            if labels[ii] == 1:
                true_positives += 1
            else:
                false_positives += 1

            ii += 1

        # Append point (1, 1)
        false_positive_rate.append(false_positives / n_negatives)
        true_positive_rate.append(true_positives / n_positives)

        return false_positive_rate, true_positive_rate

    @staticmethod
    def _get_auc(scores: np.ndarray, labels: np.ndarray) -> float:
        """ Compute the area under the ROC curve.

            Parameters
            ----------
            scores : np.ndarray of shape(n_samples,)
                The score of each sample.

            labels : np.ndarray of shape(n_samples,)
                The ground truth.

            Returns
            --------
            area : float
                The area under the curve.
        """
        scores, labels = RetrospectiveScreening._sort_scores_and_labels_for_roc(scores, labels)

        n_positives = int(np.sum(labels))
        n_negatives = labels.shape[0] - n_positives

        false_positives = 0
        true_positives = 0
        false_pos_prev = 0
        true_pos_prev = 0

        area = 0
        score_prev = -10000000

        ii = 0
        while ii < labels.shape[0]:

            if scores[ii] != score_prev:
                area += RetrospectiveScreening._trapezoid_area(false_positives, false_pos_prev,
                                                               true_positives, true_pos_prev)
                score_prev = scores[ii]
                false_pos_prev = false_positives
                true_pos_prev = true_positives

            if labels[ii] == 1:
                true_positives += 1
            else:
                false_positives += 1

            ii += 1

        area += RetrospectiveScreening._trapezoid_area(n_negatives, false_pos_prev, n_positives, true_pos_prev)
        # Scale area from n_negatives * n_positives onto the unit square
        area = area / (n_negatives * n_positives)

        assert 0 <= area <= 1

        return area

    @staticmethod
    def _trapezoid_area(x1: float, x2: float, y1: float, y2: float) -> float:
        """ Calculate the area of a trapezoid.
        """
        base = abs(x1 - x2)
        # average height
        height = abs(y1 + y2) / 2
        return base * height

    @staticmethod
    def _enrichment_data(scores: np.ndarray, labels: np.ndarray) -> Tuple[List[float], List[float]]:
        """ Get enrichment data necessary for enrichment plot and enrichment factor calculation.
        
            Parameters
            ----------
            scores : np.ndarray of shape(n_samples,)
                The score of each sample.
            
            labels : np.ndarray of shape(n_samples,)
                The ground truth.

            Returns
            --------
            screened_percentage : list of float
                The percentage of the database screened.
            
            percentage_actives_found: list of float
                The percentage of actives found.

        """
        scores, bioactivities = RetrospectiveScreening._sort_scores_and_labels_for_roc(scores, labels)

        n_molecules = bioactivities.shape[0]
        n_actives = np.sum(bioactivities)

        # Calculate % number of active molecules found in the x% of the screened database
        percentage_actives_found = []
        screened_percentage = []
        actives_counter = 0
        for ii in range(n_molecules):
            if bioactivities[ii] == 1:
                actives_counter += 1
            percentage_actives_found.append(actives_counter / n_actives)
            screened_percentage.append((ii + 1) / n_molecules)

        return screened_percentage, percentage_actives_found

    @staticmethod
    def _calculate_enrichment_factor(scores: np.ndarray, labels: np.ndarray, percentage: float) -> float:
        """ Calculate enrichment factor for the x% of the screened database 

            Parameters
            ----------
            scores : np.ndarray of shape(n_samples,)
                The score of each sample.
            
            labels : np.ndarray of shape(n_samples,)
                The ground truth.
            
            percentage : float
                Percentage of the screened database. Must be between 0 and 100
            
            Returns
            -------
            float
                The enrichment factor
        """
        screened_percentage, percentage_actives_found = RetrospectiveScreening._enrichment_data(scores, labels)
        screened_percentage = np.array(screened_percentage)
        percentage_actives_found = np.array(percentage_actives_found)

        indices_screen_per = screened_percentage <= percentage / 100
        max_enrichment_idx = np.argsort(screened_percentage[indices_screen_per])[-1]

        return percentage_actives_found[max_enrichment_idx] * 100

    def _align_molecules(self, molecules: List[Chem.Mol]) -> None:
        """ Align a list of molecules to a given pharmacophore.

        Parameters
        ----------
        molecules : list of rdkit.Chem.Mol
            List of molecules to align.

        """
        self.n_molecules += len(molecules)

        rdkit_pharmacophore, radii = self.pharmacophore.to_rdkit()
        apply_radii_to_bounds(radii, rdkit_pharmacophore)

        fdef = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
        featFactory = ChemicalFeatures.BuildFeatureFactory(fdef)

        for mol in tqdm(molecules):

            bounds_matrix = rdDistGeom.GetMoleculeBoundsMatrix(mol)
            can_match, all_matches = EmbedLib.MatchPharmacophoreToMol(mol, featFactory, rdkit_pharmacophore)
            if can_match:
                failed, _, matched_mols, _ = EmbedLib.MatchPharmacophore(all_matches,
                                                                         bounds_matrix,
                                                                         rdkit_pharmacophore,
                                                                         useDownsampling=True)
                if failed:
                    matched_mol = MolScore(0.0, mol.GetProp("_Name"), mol)
                    self.molecules.append(matched_mol)
                    continue
            else:
                matched_mol = MolScore(0.0, mol.GetProp("_Name"), mol)
                self.molecules.append(matched_mol)
                continue
            atom_match = [list(x.GetAtomIds()) for x in matched_mols]

            try:
                mol_H = Chem.AddHs(mol)
                _, embeddings, _ = EmbedLib.EmbedPharmacophore(mol_H, atom_match, rdkit_pharmacophore, count=10)
            except:
                continue

            SSDs = transform_embeddings(rdkit_pharmacophore, embeddings, atom_match)
            if len(SSDs) == 0:
                # TODO: Should we give a score of zero to molecules that don't match?
                matched_mol = MolScore(0.0, mol.GetProp("_Name"), mol)
                self.molecules.append(matched_mol)
                continue
            best_fit_index = min(enumerate(SSDs), key=itemgetter(1))[0]

            score = 1 / SSDs[best_fit_index]
            matched_mol = MolScore(score, mol.GetProp("_Name"), embeddings[best_fit_index])
            self.molecules.append(matched_mol)

    def _fingerprint_similarity(self, molecules: List[Chem.Mol]) -> None:
        """ Compute fingerprints and similarity values for a list of molecules. 

        Parameters
        ----------
        molecules : list of rdkit.Chem.Mol
            List of molecules whose similarity to the pharmacophoric fingerprint will be calculated.

        """

        self.n_molecules = len(molecules)

        for mol in tqdm(molecules):
            fingerprint = Gen2DFingerprint(mol, self._factory)
            similarity = self.similarity_fn(self.pharmacophore, fingerprint)
            mol_id = mol.GetProp("_Name")
            matched_mol = MolScore(similarity, mol_id, mol)
            self.molecules.append(matched_mol)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(n_molecules={self.n_molecules})"
