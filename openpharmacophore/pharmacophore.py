import molsysmt as msm

class Pharmacophore():

    def __init__(self):

        self.elements = set()
        self.n_elements = 0
        self.extractor = None
        self.molecular_system = None

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
