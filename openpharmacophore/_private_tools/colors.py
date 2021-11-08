# -*- coding: utf-8 -*-
"""Module with methods to work with color codes and palettes

Set of tools to work with color codes.

Notes
-----
    Most of the methods in this module are nothing but a wrapper of a method or a combination of
    methods belonging to the module `colors` of the library MatPlotLib.

.. _matplotlib colors api:
   https://matplotlib.org/stable/api/colors_api.html

"""

from matplotlib import colors as plt_colors
import numpy as np

def rgb2hex(color):

    """Method to convert a normalized RGB color code to HEX color code

    Conversion of a normalized RGB color (a list of three float numbers between 0.0 and 1.0) to HEX
    color code. This method is nothing but a wrapper of the tool `matplotlib.colors.to_hex` with the
    input argument 'keep_apha=False'.

    Parameters
    ----------
    color: :obj:`list` of :obj:`float`
        Normalized RGB color code (three float numbers between 0.0 and 1.0)

    Returns
    -------
    str
        HEX color code with the form "#RRGGBB"

    Examples
    --------

    >>> from uibcdf_stdlib.colors import rgb2hex
    >>> rgb2hex([0.0, 0.5, 0.5])
    '#008080'

    """

    return plt_colors.to_hex(color, keep_alpha=False)

def hex2rgb(color):

    """Method to convert an HEX color code to normalized RGB color code

    Conversion of an HEX color code to normalized RGB color code (a list of three float numbers
    between 0.0 and 1.0). This method is nothing but a wrapper of the tool
    `matplotlib.colors.to_rgb`.

    Parameters
    ----------
    color: str
        HEX color code

    Returns
    -------
    :obj:`list` of :obj:`float`
        Normalized RGB color code (three float numbers between 0.0 and 1.0)

    Examples
    --------

    >>> from uibcdf_stdlib.colors import hex2rgb
    >>> hex2rgb('#008080')
    (0.0, 0.5019607843137255, 0.5019607843137255)

    """

    return plt_colors.to_rgb(color)

def is_hex(color):

    """ Method to check if a color code is hexadecimal (HEX)

    This method returns the value True if the input argument corresponds to an hexadecimal color
    code.

    Parameters
    ----------
    color: :obj:
        Color code

    Returns
    -------
    bool
        True if the input color code takes the hexadecimal (HEX) form.

    Examples
    --------

    >>> from uibcdf_stdlib.colors import is_hex
    >>> is_hex('#008080')
    True
    >>> is_hex([0.0, 0.5, 0.5])
    False

    """

    output = False

    if type(color)==str:
        if (len(color)==6) or (len(color)==7 and color.startswith('#')):
            output = True

    return output

def is_rgb(color):

    """ Method to check if a color code is normalized decimal RGB.

    This method returns the value True if the input argument corresponds to a normalized decimal
    RGB color code (three float numbers in a list or tuple with values between 0.0 and 1.0).

    Parameters
    ----------
    color: :obj:
        Color code

    Returns
    -------
    bool
        True if the input color code takes the normalized decimal RGB form.

    Examples
    --------

    >>> from uibcdf_stdlib.colors import is_rgb
    >>> is_rgb('#008080')
    False
    >>> is_rgb([0.0, 0.5, 0.5])
    True

    """

    output = False

    if type(color) in [tuple, list, np.ndarray]:
        if np.shape(color)==(3,):
            if np.all([type(ii) in [float, int] for ii in color]):
                output = True

    return output

def convert(color, to_form=None):

    """ Conversion method between different formats of color codes

    This method converts color codes between different formats: normalized decimal RGB and
    hexadecimal HEX, at the moment.

    Parameters
    ----------
    color: :obj:
        Color code
    to_form: str
        Format of the output color code: 'rgb' or 'hex'.

    Returns
    -------
    :obj:
        Color code taking the format specified by the input argument `to_form`.

    Examples
    --------

    >>> from uibcdf_stdlib.colors import convert
    >>> convert('#008080', to_form='rgb')
    (0.0, 0.5019607843137255, 0.5019607843137255)
    >>> convert('#008080', to_form='hex')
    '#008080'

    """

    output = None

    if to_form=='rgb':
        if is_rgb(color):
            output = color
        elif is_hex(color):
            output = hex2rgb(color)
        else:
            raise NotImplementedError("This color format was not implemented yet")
    elif to_form=='hex':
        if is_rgb(color):
            output = rgb2hex(color)
        elif is_hex(color):
            output = color
        else:
            raise NotImplementedError("This color format was not implemented yet")

    return output

