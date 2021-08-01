from openpharmacophore import Pharmacophore

class StructuredBasedPharmacophore(Pharmacophore):

    def __init__(self, elements=[], molecular_system=None):
        super().__init__(elements=elements, molecular_system=molecular_system)
    
    
    def _from_pharmer(self, pharmacophore, load_mol_system):

        """Private method to update the attributes with those from an imported pharmer pharmacophore.

        Parameters
        ----------
        pharmacophore: :obj: str
            File or object with the pharmer pharmacophoric model.

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