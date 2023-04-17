import pyunitwizard as puw
import openpharmacophore as oph
from openpharmacophore import mol_db


def test_virtual_screening_with_a_pharmacophore():
    # We want to screen a database of molecules with a pharmacophore
    # to find potential matches. We start by defining our pharmacophore model
    pharmacophore = oph.Pharmacophore([
        oph.PharmacophoricPoint(
            feat_type="hb acceptor",
            center=puw.quantity([2.35, 1.94, -0.68], "angstroms"),
            radius=puw.quantity(1.0, "angstroms"),
        ),
        oph.PharmacophoricPoint(
            feat_type="hydrophobicity",
            center=puw.quantity([1.57, -2.53, 2.15], "angstroms"),
            radius=puw.quantity(1.0, "angstroms"),
        ),
        oph.PharmacophoricPoint(
            feat_type="aromatic ring",
            center=puw.quantity([6.19, 1.34, -1.36], "angstroms"),
            radius=puw.quantity(1.0, "angstroms"),
        ),
    ])
    # Our molecule database, in this case, is as a list of SMILES
    smiles = [
        r"[H]/N=C(\C1CCC(CC1)CNC(=O)[C@@H]2C=C(CN3N2C(=O)N(C3=O)CC(c4ccccc4)c5ccccc5)C)/N",
        r"CN[C@H](Cc1ccccc1)C(=O)N2CCC[C@H]2C(=O)NCC3CCC(CC3)N",
        r"c1ccc(cc1)S(=O)(=O)CCN2C(=O)N3CC=C[C@H](N3C2=O)C(=O)NC4CCC(CC4)c5cnc([nH]5)N",
        r"c1ccc(cc1)S(=O)(=O)CCN2C(=O)N3CC=C[C@H](N3C2=O)C(=O)NCC4CCC(CC4)N",
        r"[H]/N=C(/c1ccc(cc1)C[C@H](C(=O)N2CCCCC2)NC(=O)CNS(=O)(=O)c3ccc4ccccc4c3)\N",
        r"[H]/N=C(\c1ccc2c(c1)cc([nH]2)C(=O)N3CCC(CC3)Cc4ccccc4)/N",
        r"CCC1CCN(CC1)C(=O)[C@H](CCCNC(=[NH2+])N)NS(=O)(=O)c2cccc3c2cccc3N(C)C",
        r"C[C@H]1CN2c3c(cc(cc3NC2=S)Cl)CN1CC=C(C)C,1TVR,TB9",
        r"CC(C)c1c(n(c(n1)COC(=O)N)Cc2ccncc2)Sc3cc(cc(c3)Cl)Cl,1EP4,S11",
        r"CC(C)OC(=O)N1c2cc(ccc2NC(=S)[C@@H]1CSC)OC,1BQM,HBY",
    ]
    db = mol_db.InMemoryMolDB()
    db.from_smiles(smiles)

    # We can now create a virtual screening object to which will pass our pharmacophore
    # and our molecule db
    screener = oph.VirtualScreening([pharmacophore])
    screener.screen(db)

    assert screener.matches[0] > 0
    assert screener.fails[0] > 0
    assert all(mol.has_conformer for mol in screener.mathces[0])
    