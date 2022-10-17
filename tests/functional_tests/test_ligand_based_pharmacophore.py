import openpharmacophore as oph
import openpharmacophore.data as data


def test_ligand_based_pharmacophore_from_ligand_set():
    # We want to create a ligand-based pharmacophore for the ligands of the
    # Non-Peptide Angiotensin II Receptor

    # We start by creating a pharmacophore from a file containing the smiles
    # of the ligand
    ligand_file = data.ligands["clique_detection.smi"]
    pharmacophore = oph.load(ligand_file)
    assert isinstance(pharmacophore, oph.LigandBasedPharmacophore)

    # We confirm that the molecules were correctly loaded
    assert len(pharmacophore.ligands) == 5
    assert all([mol is not None for mol in pharmacophore.ligands])

    # We extract a pharmacophore and visualize it
    pharmacophore.extract()
    assert len(pharmacophore.pharmacophoric_points) > 0
    view = pharmacophore.show(ligands=True)
    # The view should contain a component for each ligand (5) + the pharmacophoric
    # points
    assert len(view._ngl_component_names) > 5
