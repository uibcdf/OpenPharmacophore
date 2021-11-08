from openpharmacophore import VirtualScreening
from openpharmacophore.databases import chembl, pubchem
from openpharmacophore._private_tools.exceptions import BadShapeError, OpenPharmacophoreValueError
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem

class RetrospectiveScreening(VirtualScreening):
    """ Class for performing retrospective virtual screening. 
    
        This class expects molecules classified as actives and inactives. 
        With this class pharmacophore models can be validated.

    Parameters
    ----------

    Attributes
    ----------

    """
    def __init__(self, pharmacophore, **kwargs):
        super().__init__(pharmacophore=pharmacophore, **kwargs)
        self.n_actives = 0
        self.n_inactives = 0
        self.bioactivities = None
        if self.scoring_metric == "Similarity":
            self.similarity_cutoff = 0.0

    def from_chembl_target_id(self, target_id, pIC50_threshold=6.3):
        """ Retrospective screening from bioactivity data fetched from chembl.
           
           Parameters
           ----------
           target_id : str
                ChemBl target id.
           
           pIC50_threshold : float, default=6.3
                The cuttoff value from which a molecule is considered active.
           
           """
        smiles, activity = chembl.get_training_data(target_id, pIC50_threshold)
        
        self.db = "PubChem"
        self.from_bioactivity_data(smiles, activity)

    def from_bioactivity_data(self, smiles, activity):
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
        self.n_molecules = len(smiles) 
        self.bioactivities = activity

        molecules = []
        for id, smi in smiles:
            mol = Chem.MolFromSmiles(smi)
            mol.SetProp("_Name", str(id))
            molecules.append(mol)
       
        self._screen_fn(molecules, sort=False)

    def from_pubchem_bioassay_id(self, bioassay_id):
        """ Retrospective screening from a pubchem bioassay.

            Parameters
            ----------
            bioassay_id : int
                PubChem bioassay id. 
        """
        pubchem_client = pubchem.PubChem()
        smiles, activity = pubchem_client.get_assay_training_data(bioassay_id)
        self.db = "Pubchem"
        self.from_bioactivity_data(smiles, activity)

    def AUC(self):
        """ Calculate ROC area under the curve.
        
            Returns
            -------
            area : float
                The value of the area under the curve.

            References
            ----------
            Fawcett, T. An Introduction to ROC Analysis. Pattern Recognition Letters 2006, 27, 861−874
        """
        def trapezoid_area(x1, x2, y1, y2):
            """ Calculate the area of a trapezoid.
            """
            base = abs(x1 - x2)
            # average height
            height = abs(y1 + y2) / 2
            return base * height

        if self.scoring_metric == "SSD":
            raise NotImplementedError

        scores = [x[0] for x in self.matches]
        scores = np.array(scores)

        # Sort scores in descending order and then sort labels
        indices = np.argsort(scores)[::-1]
        scores = np.sort(scores)[::-1] 
        labels = self.bioactivities[indices]

        n_positives = self.n_actives
        n_negatives = self.n_inactives

        false_positives = 0
        true_positives = 0
        false_pos_prev = 0
        true_pos_prev = 0

        area = 0
        score_prev = -10000000

        i = 0
        while i < labels.shape[0]:

            if scores[i] != score_prev:
                area += trapezoid_area(false_positives, false_pos_prev, 
                                        true_positives, true_pos_prev)
                score_prev = scores[i]
                false_pos_prev = false_positives
                true_pos_prev = true_positives
            
            if labels[i] == 1:
                true_positives += 1
            else:
                false_positives += 1

            i += 1

        area += trapezoid_area(n_negatives, false_pos_prev, n_positives, true_pos_prev)
        # Scale area from n_negatives * n_positives onto the unit square
        area = area / (n_negatives * n_positives)

        assert area >= 0 and area <= 1

        return area

    def ROC_plot(self, ax=None, label="", random_line=True):
        """ Plot the ROC curve. 
        
            Parameters
            ----------
            ax : matplotlib.axes._subplots.AxesSubplot, optional (Default = None)
                An axes object where the plot will be drawn.

            random_line : bool, default=True
                Whether to plot the line corresponding to a random classifier.

            Returns
            ----------
            ax : matplotlib.axes._subplots.AxesSubplot, optional (Default = None)
                An axes object whith the plot.

            References
            ----------
            Fawcett, T. An Introduction to ROC Analysis. Pattern Recognition Letters 2006, 27, 861−874

        """

        if self.scoring_metric == "SSD":
            raise NotImplementedError

        scores = [x[0] for x in self.matches]
        scores = np.array(scores)

        # Sort scores in descending order and then sort labels
        indices = np.argsort(scores)[::-1]
        scores = np.sort(scores)[::-1] 
        labels = self.bioactivities[indices]

        n_positives = self.n_actives
        n_negatives = self.n_inactives
        
        false_positives = 0
        true_positives = 0
        x = []
        y = []
        score_prev = -10000000

        # Calculate points for the ROC plot
        i = 0
        while i < labels.shape[0]:

            if scores[i] != score_prev:
                x.append(false_positives/n_negatives)
                y.append(true_positives/n_positives)
                score_prev = scores[i]

            if labels[i] == 1:
                true_positives += 1
            else:
                false_positives += 1

            i += 1

        # Append point (1, 1)
        x.append(false_positives/n_negatives)
        y.append(true_positives/n_positives)

        # Plot the curve

        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(x, y, label=label)
        if random_line:
            ax.plot([0, 1], [0, 1], color="black", linestyle="dashed", label="Random")
        
        ax.set_xlabel("Sensitivity")
        ax.set_ylabel("1 - Specificity")
        if label or random_line:
            ax.legend()   

        return ax

    def enrichment_factor(self, percentage):
        """ Calculate enrichment factor for the x% of the screened database 

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

        screened_percentage, percentage_actives_found = self._enrichment_data()
        screened_percentage = np.array(screened_percentage)
        percentage_actives_found = np.array(percentage_actives_found)

        indices_screen_per = screened_percentage <= percentage / 100
        max_enrichment_idx = np.argsort(screened_percentage[indices_screen_per])[-1]

        return percentage_actives_found[max_enrichment_idx] * 100

    def ideal_enrichment_factor(self, percentage):
        """ Calculate ideal enrichment factor for the x% of the screened database 

            Parameters
            ----------
            percentage : float
                Percentage of the screened database. Must be between 0 and 100
            
            Returns
            -------
            float
                The idal enrichment factor"""
    
        percentage = percentage / 100

        n_molecules = self.bioactivities.shape[0]
        ratio_actives = self.n_actives / n_molecules
        if percentage <= ratio_actives:
            return (100 / ratio_actives) * percentage
        else:
            return 100.0
    
    def enrichment_plot(self, ax=None, label="", random_line=True, ideal=False):
        """ Create an enrichment plot 
            
            Parameters
            ----------
            ax : matplotlib.axes._subplots.AxesSubplot, optional (Default = None)
                An axes object where the plot will be drawn.

            random_line : bool, default=True
                Whether to plot the line corresponding to a random classifier.

            ideal : bool, defaul=False
                Whether to plot the ideal enrichmnent curve
            
            Returns
            ----------
            ax : matplotlib.axes._subplots.AxesSubplot, optional (Default = None)
                An axes object whith the plot.
        
        """
        screened_percentage, percentage_actives_found = self._enrichment_data()
        n_molecules = self.bioactivities.shape[0]

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
    
    def _enrichment_data(self):
        """ Get enrichment data necessary for enrichment plot and enrichment factor calculation.

            Returns
            --------
            screened_percentage : list of float
                The percentage of the database screened.
            
            percentage_actives_found: list of float
                The percentage of actives found.

        """
        scores = [x[0] for x in self.matches]
        scores = np.array(scores)

        # Sort scores in descending order and then sort labels
        indices = np.argsort(scores)[::-1]
        scores = np.sort(scores)[::-1] 
        bioactivities = self.bioactivities[indices]

        n_molecules = bioactivities.shape[0]
        n_actives = self.n_actives

        # Calculate % number of active molecules found in the x% of the screened database
        percentage_actives_found = [0]
        screened_percentage = [0]
        actives_counter = 0
        for ii in range(n_molecules):
            if bioactivities[ii] == 1:
                actives_counter += 1
            percentage_actives_found.append(actives_counter / n_actives)
            screened_percentage.append(ii / n_molecules)
        
        return screened_percentage, percentage_actives_found