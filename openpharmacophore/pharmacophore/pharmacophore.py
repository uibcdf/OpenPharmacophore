from abc import ABCMeta, abstractmethod


class Pharmacophore(metaclass=ABCMeta):
    """ Abstract base class for pharmacophores.

        Parent class of LigandBasedPharmacophore and StructureBasedPharmacophore.
    """
    @abstractmethod
    def add_point(self, *args, **kwargs):
        pass

    @abstractmethod
    def remove_point(self, *args, **kwargs):
        pass

    @abstractmethod
    def remove_picked_point(self, *args, **kwargs):
        pass

    @abstractmethod
    def edit_picked_point(self, *args, **kwargs):
        pass

    @abstractmethod
    def add_point_in_picked_location(self, *args, **kwargs):
        pass

    @abstractmethod
    def extract(self, *args, **kwargs):
        pass

    @abstractmethod
    def to_json(self, *args, **kwargs):
        pass

    @abstractmethod
    def add_to_view(self, *args, **kwargs):
        pass

    @abstractmethod
    def show(self, *args, **kwargs):
        pass

    @abstractmethod
    def to_ligand_scout(self, *args, **kwargs):
        pass

    @abstractmethod
    def to_moe(self, *args, **kwargs):
        pass

    @abstractmethod
    def to_mol2(self, *args, **kwargs):
        pass

    @abstractmethod
    def to_rdkit(self, *args, **kwargs):
        pass
