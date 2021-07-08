openpharmacophore={
    'positive charge': '#5D3A00',
    'negative charge': '#F98948',
    'hb acceptor': '#D64933',
    'hb donor': '#2B303A',
    'included volume': '#0C7C59',
    'excluded volume': '#58A4B0',
    'hydrophobicity': '#BAC1B8',
    'aromatic ring': '#684E32',
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

