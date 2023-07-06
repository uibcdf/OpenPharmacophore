
def get_pharmacophore(molecular_system):

    from .load import load
    from .molecular_systems import smiles_from_pdb_id, ComplexBindingSite
    from .pharmacophore.ligand_receptor.ligand_receptor import LigandReceptorPharmacophore
    from molsysmt import convert

    protein = load(convert(molecular_system, to_form='pharm.pdb'))

    lig_ids = protein.ligand_ids()

    smiles = smiles_from_pdb_id(lig_ids[0])

    ligand = protein.get_ligand(lig_ids[0])
    ligand.fix_bond_order(smiles=smiles)
    ligand.add_hydrogens()

    protein.remove_ligand(lig_ids[0])

    protein.add_hydrogens()

    bsite = ComplexBindingSite(protein, ligand)

    pharmacophore = LigandReceptorPharmacophore(bsite, ligand)
    pharmacophore.extract()

    return protein, ligand, pharmacophore

