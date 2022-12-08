from openpharmacophore import LigandBasedPharmacophore, VirtualScreening, Pharmacophore
from openpharmacophore.io.pharmacophore_mol2 import load_mol2_pharmacophoric_points
from openpharmacophore.utils.conformers import generate_conformers
from openpharmacophore import load_from_file
import openpharmacophore.data as data
import pytest
from rdkit.Chem import MolFromSmiles


@pytest.fixture()
def ligand_based_pharmacophore():
    points = load_mol2_pharmacophoric_points(
        data.pharmacophores["elastase.mol2"]
    )
    pharma = LigandBasedPharmacophore()
    pharma.add_pharmacophore(Pharmacophore(points[0]))
    return pharma


@pytest.fixture()
def structure_based_pharmacophore():
    return load_from_file(data.pharmacophores["1M70.json"],
                          pharma_type="ligand-receptor")


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
    return ligand_based_pharmacophore.to_rdkit(0)


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
        molecules[ii] = generate_conformers(molecules[ii], 1)

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


match_mol = True


def mock_align_to_pharmacophore(*args, **kwargs):
    """ Mock align to pharmacophore that matches half the molecules"""
    global match_mol
    if match_mol:
        ssd = 2.5
        match_mol = False
        molecule = args[0]
        return molecule, ssd
    else:
        match_mol = True
        return


def test_from_list(mocker, sample_molecules, ligand_based_pharmacophore):
    mocker.patch(
        "openpharmacophore.VirtualScreening._align_to_pharmacophore",
        side_effect=mock_align_to_pharmacophore
    )
    global match_mol
    match_mol = True
    vs = VirtualScreening(ligand_based_pharmacophore)
    vs.from_list(sample_molecules, 0)

    assert len(vs.matches[0]) == 2
    assert vs.matches[0][0].mol == sample_molecules[0]
    assert vs.matches[0][0].score == 2.5
    assert vs.matches[0][1].mol == sample_molecules[2]
    assert vs.matches[0][1].score == 2.5
    assert vs.fails(0) == 1
    assert vs.num_mols(0) == 3


def assert_screening_with_files(mocker, pharmacophore, path,
                                matches, fails, n_mols, directory=False):
    mocker.patch(
        "openpharmacophore.VirtualScreening._align_to_pharmacophore",
        side_effect=mock_align_to_pharmacophore
    )
    global match_mol
    match_mol = True
    vs = VirtualScreening(pharmacophore)
    if not directory:
        vs.from_file(path, 0)
    else:
        vs.from_dir(path, 0, skip=["thrombin_ligands.sdf"])

    assert len(vs.matches[0]) == matches
    assert vs.fails(0) == fails
    assert vs.num_mols(0) == n_mols


def test_from_file_smi(mocker, ligand_based_pharmacophore):
    assert_screening_with_files(mocker, ligand_based_pharmacophore,
                                data.ligands["mols.smi"], 3, 2, 5)


def test_from_file_mol2(mocker, ligand_based_pharmacophore):
    assert_screening_with_files(mocker, ligand_based_pharmacophore,
                                data.ligands["ace.mol2"], 2, 1, 3)


def test_from_file_sdf(mocker, ligand_based_pharmacophore):
    assert_screening_with_files(mocker, ligand_based_pharmacophore,
                                data.ligands["sdf_example.sdf"], 2, 1, 3)


def test_from_dir(mocker, ligand_based_pharmacophore):
    files_dir = data.ligands["ace.mol2"]
    files_dir = "/".join(files_dir.split("/")[:-1])
    assert_screening_with_files(mocker, ligand_based_pharmacophore,
                                files_dir, 8, 8, 16, directory=True)
