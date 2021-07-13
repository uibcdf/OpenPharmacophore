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

def get_color_from_palette_for_feature(feature_name, color_palette):

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

