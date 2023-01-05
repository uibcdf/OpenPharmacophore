
class NotAPharmacophoreError(TypeError):
    """ Exception raised when an object is expected to be a pharmacophore.
    """

    def __init__(self, obj_type=None):

        self.message = ""
        if obj_type:
            self.message += f"Expected a pharmacophore. Got type {obj_type}"

        super().__init__(self.message)
