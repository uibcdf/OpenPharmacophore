from openpharmacophore import Dynophore, StructuredBasedPharmacophore, PharmacophoricPoint
import MDAnalysis
import mdtraj
import numpy as np
import pyunitwizard as puw
import pytest
from rdkit import Chem

@pytest.fixture
def mdtraj_trajectory():
    """ Loads a trajectory from an h5 file and returns an mdtraj trajectory."""
    return mdtraj.load("./openpharmacophore/data/trajectories/ERalpha.h5")

@pytest.fixture
def mdanalysis_universe():
    universe = MDAnalysis.Universe("./openpharmacophore/data/pdb/er_alpha.pdb",
                                   "./openpharmacophore/data/trajectories/er_alpha.xyz")
    return universe

def test_init_dynophore(mdtraj_trajectory, mdanalysis_universe):
    
    # Test init with mdtraj
    dynophore = Dynophore(mdtraj_trajectory)
    assert dynophore._trajectory_type == "mdt"
    assert isinstance(dynophore._trajectory, mdtraj.Trajectory)
    
    # Test init with mdAnalysis
    dynophore = Dynophore(mdanalysis_universe)
    assert dynophore._trajectory_type == "mda"
    assert isinstance(dynophore._trajectory, MDAnalysis.Universe)
    
    # Test init from file
    dynophore = Dynophore("./openpharmacophore/data/trajectories/ERalpha.h5")
    assert dynophore._trajectory_type == "mdt"
    assert isinstance(dynophore._trajectory, mdtraj.Trajectory)
    
def test_first_and_last_pharmacophore(mdtraj_trajectory):
    traj = mdtraj_trajectory
    dynophore = Dynophore(traj)
    dynophore.first_and_last_pharmacophore()
    
    assert len(dynophore.pharmacophores) == 2
    assert dynophore.n_pharmacophores == 2
    for pharmacophore in dynophore.pharmacophores:
        assert isinstance(pharmacophore, StructuredBasedPharmacophore)
        assert pharmacophore.n_pharmacophoric_points > 0
    assert dynophore.pharmacophore_indices == [0, traj.n_frames]

        
def test_pharmacophore_from_mdtraj(mdtraj_trajectory):
    dynophore = Dynophore(mdtraj_trajectory)
    
    pharmacophore_1 = dynophore._pharmacophore_from_mdtraj(0, load_ligand=True)
    pharmacophore_2 = dynophore._pharmacophore_from_mdtraj(10, load_ligand=False)
    pharmacophore_3 = dynophore._pharmacophore_from_mdtraj(15, load_mol_system=True)
    
    pharmacophores = [pharmacophore_1, pharmacophore_2, pharmacophore_3]
    for pharma in pharmacophores:
        assert isinstance(pharma, StructuredBasedPharmacophore)
        assert pharma.n_pharmacophoric_points > 0
    
    assert pharmacophore_1.ligand is not None
    assert isinstance(pharmacophore_1.ligand, Chem.Mol)
    
    assert pharmacophore_2.ligand is None
    
    assert pharmacophore_3.molecular_system is not None
    assert isinstance(pharmacophore_3.molecular_system, Chem.Mol)
    
def test_pharmacophore_from_mdanalysis():
    # TODO: Need more example trajectories to test other functions.
    # The trajectories shouldn't be heavy.
    pass

def test_pharmacophores_from_frames(mdtraj_trajectory):
    dynophore = Dynophore(mdtraj_trajectory)
    
    dynophore.pharmacophores_from_frames([10, 3, 5])
    assert dynophore.n_pharmacophores == 3
    assert dynophore.pharmacophore_indices == [10, 3, 5]
    for pharma in dynophore.pharmacophores:
        assert pharma.n_pharmacophoric_points > 0
    
    dynophore.pharmacophores_from_frames([2, 8])
    assert dynophore.n_pharmacophores == 2
    assert dynophore.pharmacophore_indices == [2, 8]
    for pharma in dynophore.pharmacophores:
        assert pharma.n_pharmacophoric_points > 0


@pytest.fixture
def pharmacophoric_points():
    """ Returns a tuple with lists of 20 pharmacophoric points of the same kind."""
    radius = puw.quantity(1.0, "angstroms")
    donor_1 = []
    for ii in range(20):
        center = np.array([1.0, 0.0, 3.0]) + ii / 10 
        donor_1.append(PharmacophoricPoint(feat_type="hb donor",
                                    center=puw.quantity(center, "angstroms"),
                                    radius=radius,
                                    atom_indices=(1,))
                       )
    donor_2 = []
    for ii in range(20):
        center = np.array([-5.0, 1.5, -2.4]) + ii / 10 
        donor_2.append(PharmacophoricPoint(feat_type="hb donor",
                                    center=puw.quantity(center, "angstroms"),
                                    radius=radius,
                                    atom_indices=(8,))
                       )
    aromatic = []
    for ii in range(20):
        center = np.array([-5.0, 1.5, -2.4]) + ii / 10 
        aromatic.append(PharmacophoricPoint(feat_type="aromatic ring",
                                    center=puw.quantity(center, "angstroms"),
                                    radius=radius,
                                    atom_indices=(2, 3, 4, 5, 6, 7))
                       )
    acceptor = []
    for ii in range(20):
        center = np.array([-0.5, 8.5, 2.4]) + ii / 10 
        acceptor.append(PharmacophoricPoint(feat_type="hb acceptor",
                                    center=puw.quantity(center, "angstroms"),
                                    radius=radius,
                                    atom_indices=(10,))
        )
    
    hydrophobic = []
    for ii in range(20):
        center = np.array([-0.5, 8.5, 2.4]) + ii / 10 
        hydrophobic.append(PharmacophoricPoint(feat_type="hydrophobicity",
                                    center=puw.quantity(center, "angstroms"),
                                    radius=radius,
                                    atom_indices=(12,))
        )
    
    return (donor_1, donor_2, aromatic, acceptor, hydrophobic)


@pytest.fixture
def pharmacophore_list(pharmacophoric_points):
    
    ligands = Chem.SDMolSupplier("./openpharmacophore/data/ligands/er_alpha_ligands.sdf")
    assert len(ligands) == 20
    
    donor_1, donor_2, aromatic, acceptor, hydrophobic = pharmacophoric_points

    # Create a list in which the first 5 pharmacophores have the same pharmacophoric_points, items
    # 8 - 12 have another 5 pharmacophores with the same pharmacophoric_points and the same for 
    # items 15 - 20 
    
    pharmacophores = []
    for ii in range(5):
        index = ii
        sb_pharmacophore = StructuredBasedPharmacophore(pharmacophoric_points=[donor_1[ii], 
                                                                    aromatic[ii], 
                                                                    acceptor[ii],
                                                                    ],
                                                    ligand=ligands[ii])
        pharmacophores.append(sb_pharmacophore)
    
    index += 1
    pharmacophores.append(StructuredBasedPharmacophore(pharmacophoric_points=[donor_1[index], 
                                                                aromatic[index], 
                                                                acceptor[index],
                                                                hydrophobic[index],
                                                                    ],
                                                    ligand=ligands[index]))
    index += 1
    pharmacophores.append(StructuredBasedPharmacophore(pharmacophoric_points=[donor_1[index],
                                                                donor_2[index], 
                                                                aromatic[index], 
                                                                acceptor[index],
                                                                hydrophobic[index],
                                                                    ],
                                                    ligand=ligands[index]))
    index += 1
    pharmacophores.append(StructuredBasedPharmacophore(pharmacophoric_points=[donor_1[index], 
                                                                aromatic[index], 
                                                                ],
                                                    ligand=ligands[index]))
    
    for ii in range(index, index + 5):
        index += 1
        sb_pharmacophore = StructuredBasedPharmacophore(pharmacophoric_points=[donor_2[ii], 
                                                                      hydrophobic[ii], 
                                                                      acceptor[ii],
                                                                      ],
                                                        ligand=ligands[ii])
        pharmacophores.append(sb_pharmacophore)
    
    index += 1
    pharmacophores.append(StructuredBasedPharmacophore(pharmacophoric_points=[donor_1[index],
                                                                donor_2[index], 
                                                                    ],
                                                    ligand=ligands[index]))
    index += 1
    pharmacophores.append(StructuredBasedPharmacophore(pharmacophoric_points=[hydrophobic[index], 
                                                                aromatic[index], 
                                                                ],
                                                    ligand=ligands[index]))
    
    for ii in range(index, index + 5):    
        sb_pharmacophore = StructuredBasedPharmacophore(pharmacophoric_points=[donor_1[ii],
                                                                  donor_2[ii], 
                                                                  hydrophobic[ii], 
                                                                  acceptor[ii],
                                                                ],
                                                        ligand=ligands[ii])    
        pharmacophores.append(sb_pharmacophore)
    
    assert len(pharmacophores) == 20
    return pharmacophores

def test_get_rpms_indices(pharmacophore_list, mdtraj_trajectory):
    
    dynophore = Dynophore(mdtraj_trajectory)
    dynophore.pharmacophores = pharmacophore_list
    dynophore.n_pharmacophores = 20
    assert len(dynophore.pharmacophores) == 20
    
    dynophore._get_unique_pharmacophoric_points(avg_coordinates=False)
    assert len(dynophore.unique_pharmacophoric_points) == 5
    
    indices = dynophore._get_rpms_indices()
    assert len(indices) == 3
    for index_list in indices:
        assert len(index_list) == 5
    
    assert indices[0] == [0, 1, 2, 3, 4]
    assert indices[1] == [8, 9, 10, 11, 12]
    assert indices[2] == [15, 16, 17, 18, 19]

def test_representative_pharmacophore_models(pharmacophore_list, mdtraj_trajectory):

    pharmacophores = pharmacophore_list

    dynophore = Dynophore(mdtraj_trajectory)
    dynophore.pharmacophores = pharmacophore_list
    dynophore.n_pharmacophores = 20

    rpms = dynophore.representative_pharmacophore_models()
    assert len(rpms) == 3
    assert rpms[0] == pharmacophores[2]
    assert rpms[1] == pharmacophores[8]
    assert rpms[2] == pharmacophores[18]
    

@pytest.fixture
def sample_dynamic_pharmacophore(mdtraj_trajectory, pharmacophoric_points):

    # Load a mock trajectory so the dynophore can be instansiated
    dynophore = Dynophore(mdtraj_trajectory)
    
    donor_1, donor_2, aromatic, acceptor, hydrophobic = pharmacophoric_points
    
    donor_1 = donor_1[0:10]
    donor_2 = donor_2[0:6]
    aromatic = aromatic[0:4]
    acceptor = acceptor[0]
    hydrophobic = hydrophobic[0]
    
    # Create 10 different pharmacophores
    pharmacophores = []
    for ii in range(10):
        if ii < len(aromatic):
            pharmacophores.append(StructuredBasedPharmacophore(
                pharmacophoric_points=[donor_1[ii], aromatic[ii]]
            ))
        else:
            pharmacophores.append(StructuredBasedPharmacophore(
                pharmacophoric_points=[donor_1[ii], donor_2[ii - 4]]
            ))
            
    pharmacophores[0].add_element(acceptor)
    pharmacophores[5].add_element(hydrophobic)
    
    dynophore.pharmacophores = pharmacophores
    dynophore.n_pharmacophores = len(pharmacophores)
    
    return dynophore

def test_get_unique_pharmacophoric_points(sample_dynamic_pharmacophore):
    dynophore = sample_dynamic_pharmacophore
    dynophore._get_unique_pharmacophoric_points(avg_coordinates=True)
    
    assert len(dynophore.unique_pharmacophoric_points) == 5
    
    # Donor 1
    donor_1 = dynophore.unique_pharmacophoric_points[0]
    center = puw.get_value(donor_1.center, "angstroms")
    assert donor_1.feature_name == "hb donor 1"
    assert donor_1.atom_indices == {1}
    assert donor_1.count == 10
    assert donor_1.frequency == 1
    assert len(donor_1.timesteps) == 10
    assert np.allclose(center, np.array([1.45, 0.45, 3.45]))
    
    # Aromatic
    aromatic = dynophore.unique_pharmacophoric_points[1]
    center = puw.get_value(aromatic.center, "angstroms")
    assert aromatic.feature_name == "aromatic ring 1"
    assert aromatic.atom_indices == {2, 3, 4, 5, 6, 7}
    assert aromatic.count == 4
    assert aromatic.frequency == 4 / 10
    assert len(aromatic.timesteps) == 4
    assert np.allclose(center, np.array([-4.85, 1.65, -2.25]))
    
    # Acceptor
    acceptor = dynophore.unique_pharmacophoric_points[2]
    center = puw.get_value(acceptor.center, "angstroms")
    assert acceptor.feature_name == "hb acceptor 1"
    assert acceptor.atom_indices == {10}
    assert acceptor.count == 1
    assert acceptor.frequency == 1 / 10
    assert np.allclose(center, np.array([-0.5, 8.5, 2.4]))

    # Donor 2
    donor_2 = dynophore.unique_pharmacophoric_points[3]
    center = puw.get_value(donor_2.center, "angstroms")
    assert donor_2.feature_name == "hb donor 2"
    assert donor_2.atom_indices == {8}
    assert donor_2.count == 6
    assert donor_2.frequency == 6 / 10
    assert np.allclose(center, np.array([-4.75, 1.75, -2.15]))
    
    # Hydrophobic
    hydrophobic = dynophore.unique_pharmacophoric_points[4]
    center = puw.get_value(hydrophobic.center, "angstroms")
    assert hydrophobic.feature_name == "hydrophobicity 1"
    assert hydrophobic.atom_indices == {12}
    assert hydrophobic.count == 1
    assert hydrophobic.frequency == 1 / 10
    assert np.allclose(center, np.array([-0.5, 8.5, 2.4]))

def test_pharmacophore_by_frequency(sample_dynamic_pharmacophore):
    dynophore = sample_dynamic_pharmacophore
    
    pharmacophore = dynophore.pharmacophore_by_frequency(0.5)
    assert pharmacophore.n_pharmacophoric_points == 2

    pharmacophoric_points = pharmacophore.pharmacophoric_points    
    assert pharmacophoric_points[0].feature_name == "hb donor 1"
    assert pharmacophoric_points[0].atom_indices == {1} 
    assert pharmacophoric_points[1].feature_name == "hb donor 2"
    assert pharmacophoric_points[1].atom_indices == {8}
    
    pharmacophore = dynophore.pharmacophore_by_frequency(0.3)
    assert pharmacophore.n_pharmacophoric_points == 3

    pharmacophoric_points = pharmacophore.pharmacophoric_points    
    assert pharmacophoric_points[1].feature_name == "aromatic ring 1"
    assert pharmacophoric_points[1].atom_indices == {2, 3, 4, 5, 6, 7} 
    
    pharmacophore = dynophore.pharmacophore_by_frequency(0.0)
    assert pharmacophore.n_pharmacophoric_points == 5
    
def test_pharmacophore_from_unique_points(sample_dynamic_pharmacophore):
    dynophore = sample_dynamic_pharmacophore
    pharmacophore = dynophore.pharmacophore_from_unique_points(["hb donor 1", "hb donor 2", "aromatic ring 1"])
    
    assert pharmacophore.n_pharmacophoric_points == 3
    assert pharmacophore.pharmacophoric_points[0].feature_name == "hb donor 1"
    assert pharmacophore.pharmacophoric_points[1].feature_name == "aromatic ring 1"
    assert pharmacophore.pharmacophoric_points[2].feature_name == "hb donor 2"

def test_pharmacophoric_point_frequency(sample_dynamic_pharmacophore):
    dataframe = sample_dynamic_pharmacophore.pharmacophoric_point_frequency()
    
    assert dataframe.shape == (5, 3)
    
    columns = list(dataframe.columns)
    assert len(columns) == 3
    assert columns[0] == "Feature Name"
    assert columns[1] == "Frequency"
    assert columns[2] == "Atoms Indices"
    
    names = dataframe["Feature Name"].tolist()
    assert names[0] == "hb donor 1"
    assert names[1] == "hb donor 2"
    assert names[2] == "aromatic ring 1"
    
    frequencies = dataframe["Frequency"].tolist()
    assert frequencies[0] == 1
    assert frequencies[1] == 6 / 10
    assert frequencies[2] == 4 / 10
    assert frequencies[3] == 1 / 10
    assert frequencies[4] == 1 / 10
    


