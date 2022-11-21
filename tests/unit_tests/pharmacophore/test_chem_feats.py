from openpharmacophore.pharmacophore.chem_feats import feature_indices


def test_feature_indices(mocker):
    mock_mol = mocker.Mock()
    mock_mol.GetSubstructMatches.side_effect = [
        ((0, 1), (2, 3)),
        ((4,),),
    ]

    patterns = ['a1aaaa1', 'a1aaaaa1']
    indices = feature_indices(patterns, mock_mol)
    assert indices == [(0, 1), (2, 3), (4,), ]
