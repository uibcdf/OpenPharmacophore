import molsysmt as msm

class Pharmacophore():

    """ Native object for pharmacophores.

    Openpharmacophore native class to store pharmacophoric models.

    Parameters
    ----------

    pharmacophore : :obj: (optional)
        File or object with pharmacophoric model. (Default: None)

    form : str (optional)
        Form of input pharmacophore: 'pharmer', 'ligandscout'. (Default: None)

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

        self.elements = elements
        self.n_elements = len(elements)
        self.molecular_system = molecular_system
        self.extractor = None

    def add_to_NGLView(self, view, color_palette='openpharmacophore'):

        """Adding the pharmacophore representation to a view (NGLWidget) from NGLView.

        Each pharmacophoric element is added to the NGLWidget as a new component.

        Parameters
        ----------
        view: :obj: `nglview.NGLWidget`
            View as NGLView widget where the representation of the pharmacophore is going to be
            added.
        color_palette: :obj: `str`, dict
            Color palette name or dictionary. (Default: 'openpharmacophore')

        Note
        ----

        Nothing is returned. The `view` object is modified in place.

        Example
        -------

        >>> import openpharmacophore as oph
        >>> import nglview as nv
        >>> pharmacophore_pharmer_file = oph.demo.pharmacophore_pharmer_file
        >>> pharmacophore = oph.Pharmacophore(pharmacophore_pharmer_file, form='pharmer')
        >>> view = nv.show_molsysmt(pharmacophore.molecular_system)
        >>> pharmacophore.add_to_NGLView(view)
        >>> view
        NGLWidget()

        """

        for element in self.elements:
            element.add_to_NGLView(view, color_palette=color_palette)

        pass

    def show(self, color_palette='openpharmacophore'):

        """Showing the pharmacophore model together with the molecular system from with it was
        extracted as a new view (NGLWidget) from NGLView.

        Parameters
        ----------
        color_palette: :obj: `str`, dict
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
        self.add_to_NGLView(view, color_palette=color_palette)

        return view

    def add_element(self, pharmacophoric_element):

        """Adding a new element to the pharmacophore.

        Parameters
        ----------
        pharmacophoric_element: :obj: `openpharmacohore.pharmacophoric_elements`
            Native object for pharmacophoric elements of any class defined in the module
            `openpharmacophore.pharmacophoric_elements`

        Returns
        -------

            The pharmacophoric element given as input argument is added to the pharmacophore
            as a new entry of the list `elements`.

        Example
        -------

        >>> import openpharmacophore as oph
        >>> pharmacophore = oph.Pharmacophore()
        >>> element = oph.pharmacophoric_elements.PositiveChargeSphere('[0,0,0] nm', '1.0 nm')
        >>> pharmacophore.add_element(element)
        >>> pharmacophore.elements
        [openpharmacophore.pharmacophoric_elements.aromatic_ring.AromaticRingSphere at 0x7fa33284a1d0>]

        """

        self.elements.append(pharmacophoric_element)
        self.n_elements +=1

        pass

    def _reset(self):

        """Private method to reset all attributes to default values.

        Note
        ----

           Nothing is returned. All attributes are set to default values.

        Example
        -------

        >>> import openpharmacophore as oph
        >>> pharmacophore = oph.Pharmacophore()
        >>> pharmacophore._from_reset()

        """

        self.elements=[]
        self.n_elements=0
        self.extractor=None
        self.molecular_system=None

    def _from_ligandscout(self, pharmacophore):

        """Private method to update the attributes with those from an imported ligandscout pharmacophore.

        Parameters
        ----------
        pharmacophore: :obj: str
            File or object with the ligandscout pharmacophoric model.

        Note
        ----

            Nothing is returned. All attributes are updated with those coming from the input pharmacophore.

        Example
        -------

        >>> import openpharmacophore as oph
        >>> pharmacophore_ligandscout_file = oph.demo.pharmacophore_ligandscout_file
        >>> pharmacophore = oph.Pharmacophore()
        >>> pharmacophore._from_pharmer(pharmacophore_ligandscout_file)

        """

        from openpharmacophore.io import from_ligandscout as _from_ligandscout
        self = _from_ligandscout(pharmacophore)

        pass

    def to_ligandscout(self, file_name=None):

        """Method to export the pharmacophore to the ligandscout compatible format.

        Parameters
        ----------
        file_name: str
            Name of file to be written with the ligandscout format of the pharmacophore.

        Note
        ----

            Nothing is returned. A new file is written.

        Example
        -------

        >>> import openpharmacophore as oph
        >>> pharmacophore = oph.demo.pharmacophore
        >>> pharmacophore.to_ligandscout('ligandscout_pharmacophore.xxx')

        """

        from openpharmacophore.io import to_ligandscout as _to_ligandscout
        return _to_ligandscout(self, file_name=file_name)

