import molsysmt as msm

class Pharmacophore():

    def __init__(self, pharmacophore=None, form=None):

        self.elements = set()
        self.n_elements = 0
        self.extractor = None
        self.molecular_system = None

        if pharmacophore is not None:

            if form=='pharmer':
                self.from_pharmer(pharmacophore)
            elif form=='ligandscout':
                self.from_ligandscout(pharmacophore)
            else:
                raise NotImplementedError

    def add_to_NGLView(self, view, color_palette='openpharmacophore'):

        for element in self.elements:
            element.add_to_NGLView(view, color_palette=color_palette)

        pass

    def show(self, color_palette='openpharmacophore'):

        view = msm.view(self.molecular_system)
        self.add_to_NGLView(view, color_palette=color_palette)

        return view

    def add_element(self, pharmacophoric_element):

        self.elements.add(pharmacophoric_element)
        self.n_elements +=1

        pass

    def from_pharmer(self, pharmacophore):

        from openpharmacophore.io import from_pharmer as _from_pharmer
        self = _from_pharmer(pharmacophore)

        pass

    def to_pharmer(self, file_name=None):

        from openpharmacophore.io import to_pharmer as _to_pharmer
        return _to_pharmer(self, file_name=file_name)

    def from_ligandscout(self, pharmacophore):

        from openpharmacophore.io import from_ligandscout as _from_ligandscout
        self = _from_ligandscout(pharmacophore)

        pass

    def to_ligandscout(self, file_name=None):

        from openpharmacophore.io import to_ligandscout as _to_ligandscout
        return _to_ligandscout(self, file_name=file_name)

