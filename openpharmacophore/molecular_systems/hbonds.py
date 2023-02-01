import mdtraj as mdt
import pyunitwizard as puw


def protein_ligand_hbonds(bsite):
    """ Get hydrogen bonds between a protein and a ligand

    Parameters
    ----------
    bsite : openPharmacophore.Protein

    Returns
    -------
    openPharmacophore.ChemFeatContainer

    """
    # TODO: obtain hydrogen bonds without using mdtraj
    traj = mdt.Trajectory(
        xyz=puw.get_value(bsite.coords),
        topology=bsite.topology.top
    )

    return mdt.baker_hubbard(traj)
