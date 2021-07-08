
class ShapeWithNoColorError(ValueError):

    def __init__(self, documentation_web=None):

        message = ('Either a feture name with a color palette or a color is needed to show'
        'this pharmacophoric shape. Check the online documentation for more information.')

        if documentation_web is not None:
            message = message[:-1] + ': {}'.format(documentation_web)

        super().__init__(message)

