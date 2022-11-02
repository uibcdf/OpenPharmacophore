import openpharmacophore as oph
import openpharmacophore.data as data
from assert_view import assert_view_contains_pharmacophore


def test_dynamic_ligand_receptor_pharmacophore():
    # We obtain pharmacophores from a md trajectory of er-alpha
    # that consists of three frames
    pharmacophore = oph.load(data.trajectories["eralpha_small.h5"])
    assert isinstance(pharmacophore, oph.LigandReceptorPharmacophore)
    # We know that the file contains a single ligand
    lig_ids = pharmacophore.receptor.ligand_ids
    assert len(lig_ids) == 1
    # The receptor already contains hydrogens
    assert pharmacophore.receptor.has_hydrogens()

    # We extract the pharmacophore
    pharmacophore.extract(lig_ids[0],
                          frames=[0, 1, 2],
                          add_hydrogens=False,
                          smiles="C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O"
                          )
    assert len(pharmacophore) == 3
    assert len(pharmacophore[0]) > 0
    assert len(pharmacophore[1]) > 0
    assert len(pharmacophore[2]) > 0

    # We inspect the ligand to see that it was correctly extracted. It should
    # have hydrogens, because the receptor has too
    assert pharmacophore.receptor.ligand.GetNumAtoms() == 44

    # We create a view of the second frame
    view = pharmacophore.show(frame=1, ligand=True, receptor=True)
    assert_view_contains_pharmacophore(view, len(pharmacophore[1]))
