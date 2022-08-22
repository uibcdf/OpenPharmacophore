from openpharmacophore.pharmacophore.color_palettes import get_color_from_palette_for_feature


def test_get_color_from_palette_for_feature():

    color = get_color_from_palette_for_feature("positive charge")
    assert color == (0.20392156862745098, 0.596078431372549, 0.8588235294117647)

    custom_palette = {
        "aromatic ring": "#32a852"
    }
    color = get_color_from_palette_for_feature("aromatic ring", custom_palette)
    assert color == (0.19607843137254902, 0.6588235294117647, 0.3215686274509804)
