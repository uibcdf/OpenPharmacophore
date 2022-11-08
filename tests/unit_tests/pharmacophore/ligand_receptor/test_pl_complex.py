import mdtraj as mdt
import numpy as np
import openpharmacophore._private_tools.exceptions as exc
import openpharmacophore.data as data
from openpharmacophore import PLComplex
import pyunitwizard as puw
import pytest
from copy import deepcopy
from rdkit import Chem
from openmm.app import PDBFile
from matplotlib.colors import to_rgb
from tst_data import estradiol


# Import for mocking
pl_module = "openpharmacophore.pharmacophore.ligand_receptor.pl_complex"
pl_class = pl_module + ".PLComplex"


@pytest.fixture()
def pl_complex():
    """ Returns a pl complex for testing. """
    pl_complex = PLComplex(data.pdb["test_with_lig.pdb"])
    pl_complex.set_ligand("EST:B")
    return pl_complex


def test_init_pl_complex(pl_complex):
    assert pl_complex.traj.n_atoms == 166
    assert pl_complex.topology.n_chains == 2


def test_is_ligand_atom_protein_atom(mocker):
    mock_atom = mocker.Mock()
    mock_atom.residue.is_water = False
    mock_atom.residue.is_protein = True
    mock_atom.residue.n_atoms = 8
    assert not PLComplex._is_ligand_atom(mock_atom)


def test_is_ligand_atom_water_atom(mocker):
    mock_atom = mocker.Mock()
    mock_atom.residue.is_water = True
    mock_atom.residue.is_protein = False
    mock_atom.residue.n_atoms = 1
    assert not PLComplex._is_ligand_atom(mock_atom)


def test_is_ligand_atom_ligand_atom(mocker):
    mock_atom = mocker.Mock()
    mock_atom.residue.is_water = False
    mock_atom.residue.is_protein = False
    mock_atom.residue.n_atoms = 15
    assert PLComplex._is_ligand_atom(mock_atom)


def test_find_ligands(pl_complex):
    assert pl_complex.ligand_ids == ["EST:B"]


def test_ligand_and_receptor_indices(pl_complex):
    pl = deepcopy(pl_complex)
    pl.ligand_and_receptor_indices()

    expected_receptor = list(range(0, 146))
    assert len(pl._receptor_indices) == len(expected_receptor)
    assert pl._receptor_indices == expected_receptor

    expected_ligand = list(range(146, 166))
    assert len(pl._lig_indices) == len(expected_ligand)
    assert pl._lig_indices == expected_ligand


def test_ligand_to_mol(pl_complex):
    pl = deepcopy(pl_complex)
    pl.ligand_and_receptor_indices()
    pl.ligand_to_mol()
    assert pl.ligand.GetNumAtoms() == 20
    assert all([
        a.GetPDBResidueInfo() is not None for a in pl.ligand.GetAtoms()
    ])


def test_ligand_to_mol_empty_indices_list(pl_complex):
    pl = deepcopy(pl_complex)
    pl.ligand_to_mol()
    assert pl.ligand.GetNumAtoms() == 20


def test_remove_ligand(pl_complex):
    pl = deepcopy(pl_complex)
    n_ligands = len(pl.ligand_ids)
    pl.ligand_and_receptor_indices()
    pl.remove_ligand()
    assert pl.traj.n_atoms == 146
    assert pl.traj.n_chains == 1
    assert pl.topology.n_atoms == 146
    assert pl.topology.n_chains == 1

    assert len(pl._lig_indices) == 0
    assert len(pl._receptor_indices) == 146
    assert len(pl.ligand_ids) == n_ligands - 1


def test_remove_ligand_empty_indices_list(pl_complex):
    pl = deepcopy(pl_complex)
    pl.remove_ligand()
    assert pl.traj.n_atoms == 146
    assert pl.topology.n_chains == 1


def test_has_hydrogens(pl_complex):
    # TODO: test with a trajectory that contains hydrogen
    assert not pl_complex.has_hydrogens()


def test_add_hydrogens(mocker, pl_complex):
    mock_model_to_traj = mocker.patch(
        pl_class + "._modeller_to_trajectory"
    )
    mock_modeller = mocker.patch(
        pl_module + ".Modeller"
    )
    pl_complex.add_hydrogens()

    mock_modeller.assert_called_once()
    _, _, kwargs = mock_modeller.mock_calls[0]
    assert len(kwargs["positions"]) == 166
    assert kwargs["topology"].getNumAtoms() == 166

    mock_model_to_traj.assert_called_once_with(mock_modeller.return_value)


@pytest.fixture()
def estradiol_mol():
    """ Returns estradiol as an rdkit.Mol object.
    """
    return estradiol


def test_fix_ligand_smiles_is_given(pl_complex, estradiol_mol):
    smiles = "CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O"
    pl = deepcopy(pl_complex)
    pl._ligand = deepcopy(estradiol_mol)
    assert not all(b.GetIsAromatic() for b in pl.ligand.GetBonds())
    assert len([a.GetSymbol() for a in pl.ligand.GetAtoms()
                if a.GetSymbol() == "H"]) == 0

    pl.fix_ligand(smiles=smiles)
    assert pl.ligand.GetNumAtoms() == 44
    assert any(b.GetIsAromatic() for b in pl.ligand.GetBonds())
    assert len([a.GetSymbol() for a in pl.ligand.GetAtoms()
                if a.GetSymbol() == "H"]) == 24


def test_fix_ligand_do_not_add_hydrogens(pl_complex, estradiol_mol):
    smiles = "CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O"
    pl = deepcopy(pl_complex)
    pl._ligand = deepcopy(estradiol_mol)
    assert not all(b.GetIsAromatic() for b in pl.ligand.GetBonds())
    assert len([a.GetSymbol() for a in pl.ligand.GetAtoms()
                if a.GetSymbol() == "H"]) == 0

    pl.fix_ligand(smiles=smiles, add_hydrogens=False)
    assert pl.ligand.GetNumAtoms() == 20
    assert any(b.GetIsAromatic() for b in pl.ligand.GetBonds())
    assert len([a.GetSymbol() for a in pl.ligand.GetAtoms()
                if a.GetSymbol() == "H"]) == 0


def test_fix_ligand_no_smiles_given(mocker, pl_complex, estradiol_mol):
    mocker.patch(
        pl_class + "._pdb_id_to_smi",
        return_value="C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O"
    )

    pl = deepcopy(pl_complex)
    pl._ligand = deepcopy(estradiol_mol)
    assert not all(b.GetIsAromatic() for b in pl.ligand.GetBonds())
    assert len([a.GetSymbol() for a in pl.ligand.GetAtoms()
                if a.GetSymbol() == "H"]) == 0

    pl.fix_ligand()
    assert pl.ligand.GetNumAtoms() == 44
    assert any(b.GetIsAromatic() for b in pl.ligand.GetBonds())
    assert len([a.GetSymbol() for a in pl.ligand.GetAtoms()
                if a.GetSymbol() == "H"]) == 24


def test_pdb_id_to_smi(mocker):
    mock_open = mocker.patch(
        pl_module + ".open",
        new=mocker.mock_open())
    mock_open.return_value.readlines.return_value = [
        "EST C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O",
        "DAO CCCCCCCCCCCC(=O)O"
    ]

    smiles = PLComplex._pdb_id_to_smi("EST:B")
    assert smiles == "C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O"
    mock_open.assert_called_once_with(data.pdb_to_smi)


def test_fix_ligand_no_smiles_found(mocker):
    mock_open = mocker.patch(
        pl_module + ".open",
        new=mocker.mock_open())
    mock_open.return_value.readlines.return_value = [
        "EST C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O",
        "DAO CCCCCCCCCCCC(=O)O"
    ]

    pl = PLComplex(data.pdb["test_with_lig.pdb"])
    pl._ligand_ids = ["ATP"]
    pl.set_ligand("ATP")
    with pytest.raises(exc.SmilesNotFoundError):
        pl.fix_ligand()


def test_fix_ligand_pdb_id_unl():
    pl = PLComplex(data.pdb["test_with_lig.pdb"])
    pl._ligand_ids = ["UNL"]
    pl.set_ligand("UNL")

    with pytest.raises(exc.SmilesNotFoundError):
        pl.fix_ligand()


def test_fix_ligand_template_and_lig_atom_number_different():
    pl = PLComplex(data.pdb["test_with_lig.pdb"])
    pl._ligand = Chem.MolFromSmiles("C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O")
    with pytest.raises(exc.DifferentNumAtomsError):
        pl.fix_ligand(smiles="CCCCCCCCCCCC(=O)O")


@pytest.fixture()
def pl_complex_no_lig():
    """ Returns a PLComplex object whose trajectory does not contain
        any ligand.

        Returns
        -------
        PLComplex
    """
    return PLComplex(data.pdb["test_no_lig.pdb"])


def test_modeller_to_trajectory():
    modeller = PDBFile(data.pdb["test_no_lig.pdb"])
    traj = PLComplex._modeller_to_trajectory(modeller)
    assert traj.n_frames == 1
    assert traj.n_atoms == 19
    assert traj.n_chains == 1
    assert traj.topology.n_atoms == 19
    assert traj.topology.n_chains == 1
    assert traj.topology.n_residues == 2

    expected_coords = np.array([
        [44.235, 80.308, 18.419],
        [43.549, 79.243, 17.706],
        [44.528, 78.252, 17.077],
        [45.699, 78.559, 16.853],
        [42.608, 79.792, 16.611],
        [43.375, 80.468, 15.608],
        [41.586, 80.762, 17.208],
        [44.030, 77.052, 16.814],
        [44.799, 76.021, 16.156],
        [44.189, 75.782, 14.791],
        [43.007, 76.033, 14.583],
        [44.721, 74.725, 16.954],
        [45.253, 74.836, 18.355],
        [46.598, 74.629, 18.621],
        [44.412, 75.147, 19.416],
        [47.098, 74.729, 19.906],
        [44.895, 75.245, 20.707],
        [46.246, 75.036, 20.946],
        [46.748, 75.126, 22.224],
    ]) / 10  # convert to nanometers

    assert np.allclose(traj.xyz, expected_coords)


def test_add_fixed_ligand(mocker, pl_complex_no_lig):
    lig_traj = mdt.load(data.pdb["estradiol.pdb"])
    mocker.patch(
        pl_module + ".mol_to_traj",
        return_value=lig_traj
    )

    pl = PLComplex(data.pdb["test_no_lig_2.pdb"])
    pl.add_fixed_ligand()

    assert pl.topology.n_atoms == 166
    assert pl._coords.shape == (1, 166, 3)
    assert len(pl._lig_indices) == 20
    assert len(pl._receptor_indices) == 146

    assert pl._ligand_ids == ["EST:B"]
    assert pl._ligand_id == "EST:B"


def test_binding_site_indices(mocker, pl_complex):
    mocker.patch(
        pl_class + ".lig_max_extent",
        return_value=puw.quantity(0.15, "nanometers")
    )
    mocker.patch(
        pl_class + ".lig_centroid",
        return_value=puw.quantity(np.array([0.0, 0.0, 0.0]), "nanometers")
    )

    pl = deepcopy(pl_complex)
    pl._coords = puw.quantity(np.array([[
        [0.3, 0.3, 0.3],
        [0.4, 0.4, 0.4],
        [0.5, 0.5, 0.5],
        [0.01, 0.01, 0.01],
        [0.02, 0.02, 0.02],
        [1., 1., 1.],
        [2., 2., 2.],
    ]]), "nanometers")
    bsite = pl.binding_site_indices(frame=0)
    assert np.all(bsite == np.array([0, 1, 2]))
    assert len(pl._lig_indices) == 20
    assert len(pl._receptor_indices) == 146


@pytest.fixture()
def expected_centroid():
    lig_coords = puw.quantity(np.array([
        [10.4106, 1.7203, 2.4775],
        [10.2995, 1.7834, 2.537],
        [10.1695, 1.7355, 2.512],
        [10.0598, 1.799, 2.5704],
        [10.1506, 1.624, 2.4274],
        [10.2621, 1.5588, 2.366],
        [10.2371, 1.4379, 2.2735],
        [10.3644, 1.3753, 2.2086],
        [10.4898, 1.3873, 2.2953],
        [10.5178, 1.5388, 2.3261],
        [10.3957, 1.6078, 2.3918],
        [10.6462, 1.5459, 2.4125],
        [10.7711, 1.4803, 2.3508],
        [10.7463, 1.3343, 2.3124],
        [10.617, 1.327, 2.2242],
        [10.6228, 1.1821, 2.1792],
        [10.7701, 1.1713, 2.1263],
        [10.8494, 1.2719, 2.2135],
        [10.961, 1.2027, 2.2746],
        [10.7379, 1.2449, 2.4419]
    ]), "nanometers")
    centroid = np.mean(lig_coords, axis=0)
    assert centroid.shape == (3,)
    return centroid


def test_lig_centroid(expected_centroid):
    pl_complex = PLComplex(data.pdb["test_with_lig.pdb"])
    pl_complex._lig_indices = list(range(146, 166))

    centroid = pl_complex.lig_centroid(frame=0)
    assert puw.is_quantity(centroid)
    assert np.allclose(expected_centroid, centroid)


def test_lig_max_extent(expected_centroid):
    pl_complex = PLComplex(data.pdb["test_with_lig.pdb"])
    pl_complex._lig_indices = list(range(146, 166))
    assert np.allclose(
        pl_complex.lig_max_extent(expected_centroid, frame=0),
        puw.quantity(0.59849, "nanometers"))


def test_ligand_feature_centroids_null_ligand_raises_error(pl_complex_no_lig):
    pl_complex = pl_complex_no_lig
    with pytest.raises(exc.NoLigandError):
        pl_complex.ligand_features("aromatic ring", 0)


def test_ligand_feature_centroids(mocker, estradiol_mol,
                                  pl_complex_no_lig):
    mock_feat_indices = mocker.patch(
        pl_class + ".feature_indices",
        return_value=[(0, 1)]
    )
    pl_complex = deepcopy(pl_complex_no_lig)
    pl_complex._coords = puw.quantity(np.array(
        [[
            [0., 0., 0.],  # Receptor
            [1., 1., 1.],
            [2., 2., 2.],  # Ligand
            [4., 4., 4.],
        ]]
    ), "nanometers")
    pl_complex._ligand = estradiol_mol
    pl_complex._lig_indices = [2, 3]

    centers, indices = pl_complex.ligand_features("hydrophobicity", frame=0)
    expected_cent = [puw.quantity(np.array([3., 3., 3.]), "nanometer")]

    assert np.all(centers[0] == expected_cent[0])
    assert indices == [[2, 3]]

    mock_feat_indices.assert_called_once()
    _, args, _ = mock_feat_indices.mock_calls[0]
    assert len(args) == 2
    assert args[0] == pl_complex.smarts_ligand["hydrophobicity"]


def test_feature_indices(mocker):
    mock_mol = mocker.Mock()
    mock_mol.GetSubstructMatches.side_effect = [
        ((0, 1), (2, 3)),
        ((4,),),
    ]

    patterns = ['a1aaaa1', 'a1aaaaa1']
    indices = PLComplex.feature_indices(patterns, mock_mol)
    assert indices == [(0, 1), (2, 3), (4,), ]


def set_receptor_feature_centroids(mocker, pl_complex_no_lig):
    mocker.patch(
        pl_class + ".feature_indices",
        return_value=[(0, 1), (2,)]
    )
    mocker.patch(
        pl_class + ".lig_centroid",
        return_value=puw.quantity(np.array([0.0, 0.0, 0.0]), "nanometers")
    )
    mocker.patch(
        pl_class + ".lig_max_extent",
        return_value=puw.quantity(0.15, "nanometers")
    )
    pl_complex = deepcopy(pl_complex_no_lig)
    return pl_complex


def test_receptor_feature_centroids(mocker, pl_complex_no_lig):
    pl_complex = set_receptor_feature_centroids(
        mocker, pl_complex_no_lig)
    pl_complex._coords = puw.quantity(np.array(
        [[
            [.2, .2, .2],
            [.4, .4, .4],
            [1., 1., 1.],
        ]]
    ), "nanometers")

    centers, indices = pl_complex.receptor_features("aromatic ring", frame=0)
    assert len(centers) == 1
    expected_cent = puw.quantity(np.array([.3, .3, .3]), "nanometer")
    assert np.allclose(centers[0], expected_cent)
    assert indices == [[0, 1]]


def test_receptor_feature_centroids_receptor_has_hydrogens(mocker, pl_complex_no_lig):
    pl_complex = set_receptor_feature_centroids(
        mocker, pl_complex_no_lig)
    pl_complex._coords = puw.quantity(np.array(
        [[
            [.2, .2, .2],
            [.1, .1, .1],
            [.3, .3, .3],
            [.4, .4, .4],
            [1., 1., 1.],
        ]]
    ), "nanometers")
    pl_complex._non_hyd_indices = [0, 3, 4]

    centers, indices = pl_complex.receptor_features("aromatic ring", frame=0)
    assert len(centers) == 1
    expected_cent = puw.quantity(np.array([.3, .3, .3]), "nanometer")
    assert np.allclose(centers[0], expected_cent)
    assert indices == [[0, 3]]


def test_receptor_features_invalid_mol_graph_raises_error(pl_complex_no_lig):
    pl_complex = deepcopy(pl_complex_no_lig)
    pl_complex._file_path = data.ligands["mols.smi"]
    with pytest.raises(exc.MolGraphError):
        pl_complex.receptor_features("aromatic ring", frame=0)


def test_hbond_indices_baker_criterion(mocker, pl_complex_no_lig):
    mock_bh = mocker.patch(
        pl_module + ".mdt.baker_hubbard",
    )
    pl_complex = pl_complex_no_lig
    pl_complex.hbond_indices(frame=0, criterion="baker")

    mock_bh.assert_called_once()
    _, args, _ = mock_bh.mock_calls[0]

    assert args[0].n_atoms == 19


@pytest.fixture()
def setup_hbonds(pl_complex_no_lig):
    pl_complex = deepcopy(pl_complex_no_lig)
    pl_complex._coords = puw.quantity(np.array([[
        [0., 0., 0.],
        [1., 1., 1.],
        [2., 2., 2.],
        [3., 3., 3.],
        [4., 4., 4.],
        [5., 5., 5.],
    ]]), "angstroms")
    pl_complex._lig_indices = [2, 3]

    hbond_indices = np.array([
        [0, 1, 2],
        [3, 4, 5]
    ])

    return pl_complex, hbond_indices


def test_hbonds_acceptors(setup_hbonds):
    pl_complex, hbond_indices = setup_hbonds

    acceptors = pl_complex.hbonds_acceptors(hbond_indices, 0)
    expected = [
        puw.quantity(np.array([2., 2., 2.]), "angstroms"),
        puw.quantity(np.array([5., 5., 5.]), "angstroms")
    ]

    assert len(acceptors) == 2
    assert np.all(acceptors[0] == expected[0])
    assert np.all(acceptors[1] == expected[1])


def test_hbonds_donors(setup_hbonds):
    pl_complex, hbond_indices = setup_hbonds

    donors = pl_complex.hbonds_donors(hbond_indices, 0)
    expected = [
        puw.quantity(np.array([0., 0., 0.]), "angstroms"),
        puw.quantity(np.array([3., 3., 3.]), "angstroms")
    ]

    assert len(donors) == 2
    assert np.all(donors[0] == expected[0])
    assert np.all(donors[1] == expected[1])


def test_interactions_view(mocker, pl_complex):
    mocker.patch(
        pl_module + ".show_mdtraj"
    )
    mock_slice = mocker.patch(
        pl_class + ".slice_traj"
    )

    hbonds = {
        "hb acceptor": [puw.quantity(np.array([1., 1., 1.]), "angstroms")]
    }
    feats = {
        "aromatic ring": [puw.quantity(np.array([2., 2., 2.]), "angstroms")]
    }
    indices = list(range(10))

    view = pl_complex.interactions_view(indices, 0, feats=[hbonds, feats])

    mock_slice.assert_called_once_with(indices, 0)
    view.add_component.assert_called_once()
    assert view.shape.add_sphere.call_count == 2
    assert view.update_representation.call_count == 2

    add_sphere_expected = [
        mocker.call([1., 1., 1.], to_rgb("#B03A2E"),
                    1.0, "hb acceptor"),
        mocker.call([2., 2., 2.], to_rgb("#F1C40F"),
                    1.0, "aromatic ring"),
    ]

    calls = view.shape.add_sphere.mock_calls
    assert calls == add_sphere_expected
    assert "ball+stick" in view.representations[0]["type"]


def test_slice_traj(pl_complex_no_lig):
    pl_complex = pl_complex_no_lig

    indices = list(range(0, 12))
    traj = pl_complex.slice_traj(indices, 0)
    assert traj.n_atoms == 7
    assert traj.n_residues == 1

    indices = list(range(6, 19))
    traj = pl_complex.slice_traj(indices, 0)
    assert traj.n_atoms == 12
    assert traj.n_residues == 1


def test_get_non_hyd_indices(pl_complex):
    pl_complex.get_non_hyd_indices()
    assert len(pl_complex._non_hyd_indices) == 166


def test_create_mol_graph_from_pdb(pl_complex_no_lig):
    pl_complex = deepcopy(pl_complex_no_lig)
    assert pl_complex._mol_graph is None

    pl_complex._create_mol_graph()
    assert pl_complex._mol_graph.GetNumAtoms() == 19


def test_create_mol_graph_from_traj_file():
    pl_complex = PLComplex(data.trajectories["pentalanine_small.gro"])
    assert pl_complex._mol_graph is None

    pl_complex._create_mol_graph()
    # Doesn't load hydrogens
    assert pl_complex._mol_graph.GetNumAtoms() == 30


def test_get_lig_conformer(estradiol_mol):
    pl = PLComplex(data.trajectories["ligand_traj.gro"])
    pl._ligand = estradiol_mol
    pl._lig_indices = list(range(20))

    mol = pl.get_lig_conformer(1)
    n_atoms = 20
    assert mol.GetNumAtoms() == n_atoms

    coords = np.zeros((n_atoms, 3))
    conf = mol.GetConformer(0)
    for ii in range(n_atoms):
        pos = conf.GetAtomPosition(ii)
        coords[ii][0] = pos.x
        coords[ii][1] = pos.y
        coords[ii][2] = pos.z

    expected = np.array([
        [30.058996, 19.800585, 23.294134],
        [28.71066, 20.192085, 23.23175],
        [27.80967, 19.451231, 22.472225],
        [26.493471, 19.81232, 22.4454],
        [28.17411, 18.286795, 21.931494],
        [29.51899, 17.838257, 21.90852],
        [29.758055, 16.657393, 21.098923],
        [31.217573, 16.271263, 20.938414],
        [32.02769, 16.571, 22.312746],
        [31.88681, 18.083992, 22.709368],
        [30.463049, 18.561434, 22.648584],
        [32.683704, 18.347963, 24.055729],
        [34.05264, 17.752693, 24.154547],
        [34.150806, 16.258833, 23.64093],
        [33.56477, 16.15493, 22.23642],
        [33.983757, 14.81229, 21.751375],
        [35.410973, 14.556562, 22.259155],
        [35.55814, 15.64381, 23.291952],
        [36.178825, 15.285156, 24.4787],
        [33.49326, 15.415669, 24.809036]])
    assert np.allclose(coords, expected, rtol=0.1)
