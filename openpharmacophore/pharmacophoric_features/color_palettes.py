openpharmacophore={
    'positive_charge': '#5D3A00',
    'negative_charge': '#F98948',
    'hb_acceptor': '#D64933',
    'hb_donor': '#2B303A',
    'inclusion_volume': '#0C7C59',
    'exclusion_volume': '#58A4B0',
    'hydrophobic': '#BAC1B8',
    'aromatic': '#684E32',
}

def get_color_from_palette_for_feature(feature, color_palette):

    if type(color_palette)==str:
        try:
            color_palette = globals()[color_palette]
        except:
            raise InputArgumentError('color_palette')

    try:
        color = color_palette[feature]
    except:
        raise InputArgumentError('feature')

    return color

