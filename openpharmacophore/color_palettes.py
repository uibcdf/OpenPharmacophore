"""Module with objects and methods to choose and define the color code to represent pharmacophoric
features when a pharmacophore is shown.

This module contains predefined color palettes to represent pharmacophoric features, as well as a
method to extract color codes from these palettes or from a user defined one as a dictionary.

Attributes
----------
openpharmacophore: dict
    Native color palette.

Todo
----
Some color palettes can be included here. The native 'openpharmacophore' palette should be
redefined: colors are temporary.

"""

openpharmacophore={
    'positive charge': '#E1B07E',
    'negative charge': '#A5F8D3',
    'hb acceptor': '#F13030',
    'hb donor': '#5B618A',
    'included volume': '#109648',
    'excluded volume': '#14110F',
    'hydrophobicity': '#9EADC8',
    'aromatic ring': '#D6D84F',
}

def get_color_from_palette_for_feature(feature_name, color_palette='openpharmacophore'):

    """Getting the color code of a pharmacophoric feature from a color palette.

    A color palette is a Python dictionary where keys are the farmacophoric feature names and
    the corresponding values are hexadecimal or normalized RGB color codes. For example:

        my_color_palette={
            'positive charge': '#E1B07E',
            'negative charge': '#A5F8D3',
            'hb acceptor': '#F13030',
            'hb donor': '#5B618A',
            'included volume': '#109648',
            'excluded volume': '#14110F',
            'hydrophobicity': '#9EADC8',
            'aromatic ring': '#D6D84F',
        }

    Some colorpalettes are defined already in the module
    `openpharmacophore.pharmacophoric_elements.features.color_palettes`.

    Parameters
    ----------
    feature_name: str
        Feature name: 'positive charge', 'negative charge', 'hb acceptor', 'hb donor',
        'included volume', 'excluded volume', 'hydrophobicity' or 'aromatic ring'.
    color_palette: :obj: `str`, dict
        Dictionary or color palette name predefined already in the module `openpharmacophore.pharmacophoric_elements.features.color_palettes`. (Default: 'openpharmacophore')

    Examples
    -------
    >>> import openpharmacophore as oph
    >>> my_color_palette={
    ... 'positive charge': '#E1B07E',
    ... 'negative_charge': '#A5F8D3',
    ... 'hb acceptor': '#F13030',
    ... 'hb donor': '#5B618A',
    ... 'included volume': '#109648',
    ... 'excluded volume': '#14110F',
    ... 'hydrophobicity': '#9EADC8',
    ... 'aromatic ring': '#D6D84F',
    ... }
    >>> oph.pharmacophoric_elements.features.color_palettes.get_color_from_palette_for_feature('hb donor', my_color_palette)
    '#5B618A'

    >>> import openpharmacophore as oph
    >>> oph.pharmacophoric_elements.features.color_palettes.get_color_from_palette_for_feature('positive charge', 'openpharmacophore')
    '#E1B07E'

    Returns
    -------
    :obj: str, :obj:`list` of :obj:`float`
        Color code for input feature in color palette.

    """

    if type(color_palette)==str:
        try:
            color_palette = globals()[color_palette]
        except:
            raise InputArgumentError('color_palette')

    try:
        color = color_palette[feature_name]
    except:
        raise InputArgumentError('feature_name')

    return color

