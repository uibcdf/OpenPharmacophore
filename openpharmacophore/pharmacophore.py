from openpharmacophore._private_tools.exceptions import InvalidFeatureError
import pyunitwizard as puw
import nglview as nv
from openpharmacophore.pharmacophoric_elements.features.color_palettes import get_color_from_palette_for_feature

class Pharmacophore():

    """ Native object for pharmacophores.

    Parent class of LigandBasedPharmacophore, StrucutredBasedPharmacophore and Dynophore

    Openpharmacophore native class to store pharmacophoric models.

    Parameters
    ----------

    elements : :obj:`list` of :obj:`openpharmacophore.pharmacoforic_elements`
        List of pharmacophoric elements

    molecular_system : an rdkit.Chem.rdchem.Mol or a list of rdkit.Chem.rdchem.Mol
        Set of ligands from which this pharmacophore was extracted.

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

    def add_to_NGLView(self, view, palette='openpharmacophore'):

        """Adding the pharmacophore representation to a view (NGLWidget) from NGLView.

        Each pharmacophoric element is added to the NGLWidget as a new component.

        Parameters
        ----------
        view: :obj: `nglview.NGLWidget`
            View as NGLView widget where the representation of the pharmacophore is going to be
            added.
        palette: :obj: `str`, dict
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

        # TODO: Openpharmacophore palette is not working

        for i, element in enumerate(self.elements):
            # Add Spheres
            center = puw.get_value(element.center, to_unit="angstroms").tolist()
            radius = puw.get_value(element.radius, to_unit="angstroms")
            # feature_color = get_color_from_palette_for_feature(element.feature_name, color_palette=palette)
            feature_color = feature_colors[element.feature_name]
            label = f"{element.feature_name}_{i}"
            view.shape.add_sphere(center, feature_color, radius, label)
            # Add vectors
            if element.has_direction:
                end_arrow = puw.get_value(element.center + 2 * radius * puw.quantity(element.direction, "angstroms"), to_unit='angstroms').tolist()
                label = f"{element.feature_name}_vector"
                view.shape.add_arrow(center, end_arrow, feature_color, 0.2, label)



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

        view = nv.NGLWidget()
        self.add_to_NGLView(view, palette=palette)

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
    
    def remove_element(self, element_indx):

        """Remove an element from the pharmacophore.

        Parameters
        ----------
        element_inx: int
            Index of the element to be removed

        Returns
        -------
            The pharmacophoric element given as input argument is removed from the pharmacophore.
        """

        self.elements.pop(element_indx)
        self.n_elements -=1

    def remove_feature(self, feat_type):

        """Remove an especific feature type from the pharmacophore elements list

        Parameters
        ----------
        feat_type: str
            Name or type of the feature to be removed

        Returns
        ------
            The pharmacophoric elements of the feature type given as input argument 
            are removed from the pharmacophore.
        """
        feats = [
            "aromatic ring",
            "hydrophobicity",
            "hb acceptor",
            "hb donor",
            "included volume",
            "excluded volume",
            "positive charge",
            "negative charge",
        ]    
        if feat_type not in feats:
            raise InvalidFeatureError(f"Cannot remove feature. \"{feat_type}\" is not a valid feature type")
        temp_elements = [element for element in self.elements if element.feature_name != feat_type]
        if len(temp_elements) == self.n_elements: # No element was removed
            raise InvalidFeatureError(f"Cannot remove feature. The pharmacophore does not contain any {feat_type}")
        self.elements = temp_elements
        self.n_elements = len(self.elements)

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
        # TODO: complete function

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
        # TODO: complete function

        from openpharmacophore.io import to_ligandscout as _to_ligandscout
        return _to_ligandscout(self, file_name=file_name)

