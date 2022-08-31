from openpharmacophore.pharmacophore.pharmacophore import Pharmacophore


class LigandBasedPharmacophore(Pharmacophore):

    def __init__(self, ligands):
        self._points = []

    @property
    def pharmacophoric_points(self):
        return self._points

    @property
    def num_points(self, *args, **kwargs):
        return len(self._points)

    @classmethod
    def from_file(cls, file_name):
        pass

    def add_point(self, *args, **kwargs):
        pass

    def remove_point(self, *args, **kwargs):
        pass

    def remove_picked_point(self, *args, **kwargs):
        pass

    def edit_picked_point(self, *args, **kwargs):
        pass

    def add_point_in_picked_location(self, *args, **kwargs):
        pass

    def to_json(self, *args, **kwargs):
        pass

    def add_to_view(self, *args, **kwargs):
        pass

    def show(self, *args, **kwargs):
        pass
