from openpharmacophore.pharmacophoric_elements.positive_charge import PositiveChargeGaussianKernel
from openpharmacophore import Pharmacophore
import molsysmt as msm

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
    
    def from_pdb():
        # TODO: complete function
        pass

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

        Example
        -------

        >>> import openpharmacophore as oph
        >>> import nglview as nv
        >>> pharmacophore_pharmer_file = oph.demo.pharmacophore_pharmer_file
        >>> pharmacophore = oph.Pharmacophore(pharmacophore_pharmer_file, form='pharmer')
        >>> view = pharmacophore.show()
        >>> view
        NGLWidget()

        """

        view = msm.view(self.molecular_system, standardize=False)
        self.add_to_NGLView(view, palette=palette)
        return view
    
