import openpharmacophore as oph
import openpharmacophore.data as data
import pytest
from rdkit import Chem


def has_hydrogens(mol):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "H":
            return True
    return False


def test_ligand_based_pharmacophore_ligands_preparation():
    # We want to prepare the ligands of Thrombin for pharmacophore extraction

    # We start by creating a LigandBasedPharmacophore object from a list of smiles representing
    # Thrombin ligands
    ligands_smi = [
        r"[H]/N=C(\C1CCC(CC1)CNC(=O)[C@@H]2C=C(CN3N2C(=O)N(C3=O)CC(c4ccccc4)c5ccccc5)C)/N",
        r"CN[C@H](Cc1ccccc1)C(=O)N2CCC[C@H]2C(=O)NCC3CCC(CC3)N",
        r"c1ccc(cc1)S(=O)(=O)CCN2C(=O)N3CC=C[C@H](N3C2=O)C(=O)NC4CCC(CC4)c5cnc([nH]5)N",
        r"c1ccc(cc1)S(=O)(=O)CCN2C(=O)N3CC=C[C@H](N3C2=O)C(=O)NCC4CCC(CC4)N",
    ]
    pharmacophore = oph.LigandBasedPharmacophore()
    pharmacophore.load_ligands_from_smi(ligands_smi)

    # We confirm that the molecules were correctly loaded
    assert len(pharmacophore.ligands) == 4
    assert all([mol is not None for mol in pharmacophore.ligands])

    # We proceed to add hydrogens to the ligands
    pharmacophore.add_hydrogens(ligands="all")
    assert all([has_hydrogens(mol) for mol in pharmacophore.ligands])

    # We generate conformers for each
    pharmacophore.generate_conformers(ligands="all", n_confs=1)
    assert all([mol.GetNumConformers() == 1 for mol in pharmacophore.ligands])

    # Now we want to find chemical features in the ligands
    pharmacophore.find_chem_feats()
    # All ligands have aromatic rings, hydrophobic areas, acceptors and donors
    assert len(pharmacophore.feats) == 4
    assert all(["R" in f for f in pharmacophore.feats])
    assert all(["A" in f for f in pharmacophore.feats])
    assert all(["D" in f for f in pharmacophore.feats])
    assert all(["H" in f for f in pharmacophore.feats])
    # Finally we create a 2D representation of the ligands with their
    # chemical features highlighted
    pharmacophore.draw((300, 280))


def read_sdf(file_path):
    """ Reads a sdf file that can contain multiple conformers
        for a molecule.

        Parameters
        ----------
        file_path : str

        Returns
        -------
        list : [rdkit.Chem.Mol]
    """
    supp = Chem.SDMolSupplier(file_path, removeHs=False)
    molecules = {}

    for mol in supp:
        name = mol.GetProp("_Name")
        try:
            molecules[name].AddConformer(mol.GetConformer(), assignId=True)
        except KeyError:
            molecules[name] = mol

    return list(molecules.values())


def test_ligand_based_pharmacophore_extraction():
    # We want to extract a ligand based pharmacophore for thrombin.
    # We start by loading the ligands from a sdf file
    pharmacophore = oph.LigandBasedPharmacophore()
    pharmacophore.ligands = read_sdf(data.ligands["thrombin_ligands.sdf"])
    assert len(pharmacophore.ligands) == 7
    # We extract pharmacophores of 3 points and visualize them
    pharmacophore.find_chem_feats()
    pharmacophore.extract(
        n_points=5, min_actives=5, max_pharmacophores=10
    )
    assert len(pharmacophore) > 0

    # We inspect the features of a pharmacophore. We expect thrombin pharmacophore
    # to have an aromatic ring, at leas one hydrophobic and a positive charge
    feat_names = [p.feature_name for p in pharmacophore[0]]
    assert "hydrophobicity" in feat_names
    assert "aromatic ring" in feat_names
    assert "positive charge" in feat_names

    # We visualize the best scoring pharmacophore
    view = pharmacophore.show(index=0, ligands=True)
    # The view should contain a component for each ligand (5) + the pharmacophoric
    # points
    assert len(view._ngl_component_names) > 5
