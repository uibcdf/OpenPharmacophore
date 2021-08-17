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
    
    
    def _from_pharmer(self, pharmacophore, load_mol_system=False):

        """Private method to update the attributes with those from an imported pharmer pharmacophore.

        Parameters
        ----------
        pharmacophore: :obj: str
            File or object with the pharmer pharmacophoric model.

        load_mol_sys: bool
            If true loads the molecular system associated to the pharmacophore. (Default: False)

        Note
        ----

            Nothing is returned. All attributes are updated with those coming from the input pharmacophore.

        Example
        -------

        >>> import openpharmacophore as oph
        >>> pharmacophore_pharmer_file = oph.demo.pharmacophore_pharmer_file
        >>> pharmacophore = oph.Pharmacophore()
        >>> pharmacophore._from_pharmer(pharmacophore_pharmer_file)

        """

        from openpharmacophore.io import from_pharmer as _from_pharmer
        tmp_pharmacophore = _from_pharmer(pharmacophore, load_mol_system)
        self._reset()
        self.elements = tmp_pharmacophore.elements
        self.n_elements = tmp_pharmacophore.n_elements
        self.molecular_system = tmp_pharmacophore.molecular_system

        pass

    def to_pharmer(self, file_name=None):

        """Method to export the pharmacophore to the pharmer compatible format.

        Parameters
        ----------
        file_name: str
            Name of file to be written with the pharmer format of the pharmacophore.

        Note
        ----

            Nothing is returned. A new file is written.

        Example
        -------

        >>> import openpharmacophore as oph
        >>> pharmacophore = oph.demo.pharmacophore
        >>> pharmacophore.to_pharmer('pharmer_pharmacophore.json')

        """

        from openpharmacophore.io import to_pharmer as _to_pharmer
        return _to_pharmer(self, file_name=file_name)
