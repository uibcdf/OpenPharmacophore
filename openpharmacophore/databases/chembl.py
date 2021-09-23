from chembl_webresource_client.new_client import new_client
import numpy as np
import pandas as pd
from tqdm.auto import tqdm
import os

def get_ro5_dataset( download_dir):
    """Download subset of molecules that do not violate Lipinky's rule of five.
    Molecules are stored in smi files with its smiles and ChemblId
        
        Parameters
        ----------
        download_dir: str
            Directory where the files will be saved.
        
        Notes
        ---------
        Nothing is returned. New files are written.

    """
    molecules_api = new_client.molecule
    no_violations = molecules_api.filter(molecule_properties__num_ro5_violations=0).only(
        "molecule_chembl_id",
        "molecule_structures"
    )

    # Number of molecules per file
    n_molecules = 10000
    print("Downloading molecules...")
    
    file_number = 1
    file_name =  f"mols_{file_number:02d}.smi"
    file = open(os.path.join(download_dir, file_name), "w")

    counter = 0
    for mol in tqdm(no_violations):

        # Create a new file when n_molecules have been written
        if counter % n_molecules == 0 and counter != 0:
            file.close()
            file_number += 1
            file_name =  f"mols_{file_number:02d}.smi"
            file = open(os.path.join(download_dir, file_name), "w")
            
        if mol["molecule_structures"] is not None:
            counter += 1
            mol_id = mol["molecule_chembl_id"]
            smiles = mol["molecule_structures"]["canonical_smiles"]
            line = smiles + " " + mol_id + "\n"
            file.write(line)

    file.close()


def get_bioactivity_data(target_chembl_id):
    """Get bioactivity data for a given target.
    
        Parameters
        ----------
        target_chembl_id: str
            The target chembl id.
        
        Returns
        -------
        bioactivities_df: pandas.DataFrame
            A dataframe with the following columns: ChemblID, Smiles,
            and pIC50.
    """
    activity_api = new_client.activity
    compounds_api = new_client.molecule

    bioactivities = activity_api.filter(
        target_chembl_id=target_chembl_id,
        type="IC50",
        relation="=",
        assay_type="B"
    ).only(
        "molecule_chembl_id",
        "standard_units",
        "standard_value",
    )

    n_records = len(bioactivities)

    mol_chembl_id = np.empty(n_records, dtype=object)
    IC50 = np.zeros(n_records)
    units = np.empty(n_records, dtype=object)

    print("Downloading bioactivity data...")
    for i, bioact in enumerate(tqdm(bioactivities)):
        mol_chembl_id[i] = bioact["molecule_chembl_id"]
        IC50[i] = bioact["standard_value"]
        units[i] = bioact["standard_units"]

    bioactivities_dict = {
        "ChemblID": mol_chembl_id,
        "IC50": IC50,
        "Units": units 
    }

    bioactivities_df = pd.DataFrame().from_dict(bioactivities_dict)

    # Clean bioactivities_df
    bioactivities_df = bioactivities_df.dropna(axis=0, how="any")
    bioactivities_df = bioactivities_df[bioactivities_df["Units"] == "nM"]
    # Remove duplicate elements. Keep the mean of the IC50
    bioactivities_df = bioactivities_df.groupby("ChemblID").mean().reset_index()

    # Get compounds smiles
    compounds = compounds_api.filter(
            molecule_chembl_id__in=list(bioactivities_df["ChemblID"])
        ).only("molecule_chembl_id", 
            "molecule_structures")

    comp_id = []
    smiles = []
    print("Downloading compounds smiles...")
    for comp in tqdm(compounds):
        if comp["molecule_structures"] is not None:
            comp_id.append(comp["molecule_chembl_id"])
            smiles.append(comp["molecule_structures"]["canonical_smiles"])

    compounds_dict = {
        "ChemblID": comp_id,
        "Smiles": smiles
    }

    compounds_df = pd.DataFrame.from_dict(compounds_dict)
    bioactivities_df = pd.merge(
                bioactivities_df,
                compounds_df,
                on="ChemblID",
            )
    bioactivities_df["pIC50"] = bioactivities_df.apply(
        lambda x: 9 - np.log10(x.IC50), 
        axis=1)
    bioactivities_df = bioactivities_df.drop("IC50", axis=1)

    return bioactivities_df
    
def get_training_data( target_chembl_id, pIC50_threshold=6.3):
    """Get a list of active and inactive compounds for a given target.

        Parameters
        ----------
        target_chembl_id: str
            The target chembl id.

        pIC50_threshold: float
            The cuttoff value from which a molecule is considered active.
        
        Returns
        -------
        actives: 2-tuple
            First element of the tuple is a list of Chembl ids. Second element
            is a list of smiles for the active compounds
        
        inactives: 2-tuple
            First element of the tuple is a list of Chembl ids. Second element
            is a list of smiles for the inactive compounds
    """
    bioactivities_df = get_bioactivity_data(target_chembl_id=target_chembl_id)

    actives_df = bioactivities_df[bioactivities_df["pIC50"] >= pIC50_threshold]
    inactives_df = bioactivities_df[bioactivities_df["pIC50"] < pIC50_threshold]

    actives = (actives_df["ChemblID"].tolist(), actives_df["Smiles"].tolist())
    inactives = (inactives_df["ChemblID"].tolist(), inactives_df["Smiles"].tolist())

    return actives, inactives
