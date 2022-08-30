import openpharmacophore.data as data
from openpharmacophore import LigandBasedPharmacophore, PharmacophoricPoint
from openpharmacophore.utils.conformers import generate_conformers
import pytest
import pyunitwizard as puw
from rdkit import Chem


@pytest.fixture
def ligand_based_pharmacophore():
    """ Returns a pharmacophore model for the
        Angiotensin II Receptor antagonists"""
    smiles = Chem.SmilesMolSupplier(data.ligands["mols"], titleLine=0)
    ligands = [mol for mol in smiles]
    assert len(ligands) == 5
    radius = puw.quantity(1.0, "angstroms")
    points = [
        PharmacophoricPoint("hb acceptor",
                            puw.quantity([-3.2099, 0.6359, 2.4323], "angstroms"),
                            radius),
        PharmacophoricPoint("hb acceptor",
                            puw.quantity([-1.5734, 4.2674, 0.3026], "angstroms"),
                            radius),
        PharmacophoricPoint("aromatic ring",
                            puw.quantity([1.0191, 0.7249, -0.9843], "angstroms"),
                            radius),
        PharmacophoricPoint("aromatic ring",
                            puw.quantity([-2.6875, 1.3437, 1.7495], "angstroms"),
                            radius),
        PharmacophoricPoint("hb donor",
                            puw.quantity([4.146, -1.3919, 1.8679], "angstroms"),
                            radius),
        PharmacophoricPoint("hydrophobicity",
                            puw.quantity([-4.8432, -1.7925, -0.4132], "angstroms"),
                            radius),
    ]

    return LigandBasedPharmacophore(points, ligands)


@pytest.mark.skip(reason="Extraction of ligand based pharmacophores has not been implemented yet")
def test_from_ligand_list():
    assert False, "Complete me!"


@pytest.mark.skip(reason="Extraction of ligand based pharmacophores has not been implemented yet")
def test_from_ligand_file():
    assert False, "Complete me!"


def test_single_ligand_pharmacophore():
    phenol = Chem.MolFromSmiles("C1=CC=C(C=C1)O")
    phenol = generate_conformers(phenol, 1)
    pharmacophore = LigandBasedPharmacophore.single_ligand(phenol,
                                                           features=["hb donor", "aromatic ring"])
    assert len(pharmacophore) == 2
    assert pharmacophore[0].feature_name == "hb donor"
    assert pharmacophore[1].feature_name == "aromatic ring"


def test_show(ligand_based_pharmacophore):

    view_no_ligands = ligand_based_pharmacophore.show(show_ligands=False)
    assert len(view_no_ligands._ngl_component_ids) == 6
    for component in view_no_ligands._ngl_component_names:
        assert component == "nglview.shape.Shape"

    view_with_ligands = ligand_based_pharmacophore.show()
    # 6 pharmacophoric points + 5 ligands
    assert len(view_with_ligands._ngl_component_ids) == 11
    for ii in range(5):
        assert view_with_ligands._ngl_component_names[ii] == "nglview.adaptor.RdkitStructure"
    for ii in range(5, 11):
        assert view_with_ligands._ngl_component_names[ii] == "nglview.shape.Shape"
