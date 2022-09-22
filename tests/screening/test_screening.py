from openpharmacophore import LigandBasedPharmacophore, StructureBasedPharmacophore, VirtualScreening
from openpharmacophore.utils.conformers import generate_conformers
import openpharmacophore.data as data
import pytest
from rdkit.Chem import MolFromSmiles


@pytest.fixture()
def ligand_based_pharmacophore():
    pharmacophore = LigandBasedPharmacophore([])
    pharmacophore.from_file(data.pharmacophores["elastase"])
    return pharmacophore


@pytest.fixture()
def structure_based_pharmacophore():
    pharmacophore = StructureBasedPharmacophore(None)
    pharmacophore.from_file(data.pharmacophores["1M70"])
    return pharmacophore


def test_init_virtual_screening_with_lbp(ligand_based_pharmacophore):
    vs = VirtualScreening(ligand_based_pharmacophore)
    assert len(vs.pharmacophores) == 1
    assert len(vs.pharmacophores[0].getFeatures()) == 4
    assert len(vs.matches) == 1


def test_init_virtual_screening_with_sbp(structure_based_pharmacophore):
    vs = VirtualScreening(structure_based_pharmacophore)
    assert len(vs.pharmacophores) == 1
    assert len(vs.pharmacophores[0].getFeatures()) == 5
    assert len(vs.matches) == 1


@pytest.fixture()
def rdkit_pharmacophore(ligand_based_pharmacophore):
    return ligand_based_pharmacophore.to_rdkit()


def test_align_to_pharmacophore_mol_features_dont_match(mocker, rdkit_pharmacophore):
    mock_embed = mocker.patch(
        "openpharmacophore.screening.screening.EmbedLib"
    )
    mock_embed.MatchPharmacophoreToMol.return_value = (False, [])
    mol = MolFromSmiles("C")
    mol_and_score = VirtualScreening._align_to_pharmacophore(mol, rdkit_pharmacophore)
    assert mol_and_score is None


def test_align_to_pharmacophore_mol_doesnt_match(mocker, rdkit_pharmacophore):
    mock_embed = mocker.patch(
        "openpharmacophore.screening.screening.EmbedLib"
    )
    mock_embed.MatchPharmacophoreToMol.return_value = (True, [])
    mock_embed.MatchPharmacophore.return_value = (True, None, None, None)

    mol = MolFromSmiles("C")
    mol_and_score = VirtualScreening._align_to_pharmacophore(mol, rdkit_pharmacophore)
    assert mol_and_score is None


@pytest.fixture()
def sample_molecules():
    """ Returns a list of three molecules each with one conformer.
    """
    molecules = [
        MolFromSmiles("c1ccccc1"),
        MolFromSmiles("n1ccccc1"),
        MolFromSmiles("o1cccc1")
    ]
    for ii in range(len(molecules)):
        molecules[ii] = generate_conformers(molecules[ii])

    return molecules


def test_align_to_pharmacophore_matching_mol(mocker, rdkit_pharmacophore,
                                             sample_molecules):
    mock_feature = mocker.Mock()
    mock_feature.GetAtomIds.return_value = (1, 2)
    matched_features = [mock_feature] * 2
    mock_embed = mocker.patch(
        "openpharmacophore.screening.screening.EmbedLib")
    mock_embed.MatchPharmacophoreToMol.return_value = (True, [])
    mock_embed.MatchPharmacophore.return_value = (False, None, matched_features, None)
    mock_embed.EmbedPharmacophore.return_value = (None, sample_molecules, 0)

    mock_transform = mocker.patch(
        "openpharmacophore.VirtualScreening._transform_embeddings")
    mock_transform.return_value = [1.5, 1.0, 2.0]

    mol = MolFromSmiles("C")
    assert VirtualScreening._align_to_pharmacophore(mol, rdkit_pharmacophore) == (
        sample_molecules[1], 1.0)


def test_transform_embeddings(mocker, rdkit_pharmacophore):
    mock_trans_matrix = mocker.patch(
        "openpharmacophore.VirtualScreening._transform_matrix"
    )
    mocker.patch("rdkit.Chem.rdMolTransforms.TransformConformer")
    mock_trans_matrix.return_value = (1.5, None)

    molecule = MolFromSmiles("c1ccccc1")
    embeddings = [
        generate_conformers(molecule, 1),
        generate_conformers(molecule, 1)
    ]
    assert VirtualScreening._transform_embeddings(
        rdkit_pharmacophore, embeddings, [[1, 2], [2, 3]]
    ) == [1.5, 1.5]


def test_transform_matrix(mocker, sample_molecules):
    mock_alignment = mocker.patch(
        "rdkit.Numerics.rdAlignment.GetAlignmentTransform")
    # We are not going to test the return value. Just that it
    # doesn't raise an error
    mock_alignment.return_value = (1.0, None)

    conformer = sample_molecules[0].GetConformer()
    atom_match = [[1, 2], [0, 3]]
    assert VirtualScreening._transform_matrix(
        None, conformer, atom_match
    ) == (1.0, None)
