from openpharmacophore import Ligand, Pharmacophore
from openpharmacophore.pharmacophore.ligand_based.common_pharmacophore import CommonPharmacophoreFinder


class LigandBasedPharmacophore:
    """ Class to store and extract pharmacophores from a set of ligands.

        Parameters
        ----------
        ligands : list[Ligand]

    """
    def __init__(self, ligands):
        self._ligands = ligands
        self._pharmacophores = []  # type: list[Pharmacophore]

    def extract(self, n_points, min_actives=None, max_pharmacophores=None, *args, **kwargs):
        """ Finds and scores common pharmacophores from a set of ligands.

           Parameters
           ----------
           n_points : int
               Extracted pharmacophores will have this number of pharmacophoric
               points.

           min_actives : int, optional
               Number of ligands that must match a common pharmacophore.

           max_pharmacophores : int, optional
               Maximum number of pharmacophores to return. If set to null
               all found pharmacophores will be returned.
        """
        if min_actives is None:
            min_actives = len(self._ligands)

        chem_feats = []
        for lig in self._ligands:
            ligand_feats = []
            for conf in range(lig.n_conformers):
                ligand_feats.append(lig.get_chem_feats(conf, indices=True))
            chem_feats.append(ligand_feats)

        extractor = CommonPharmacophoreFinder(n_points, min_actives, max_pharmacophores, *args, **kwargs)
        return extractor(chem_feats)

    def __len__(self):
        return len(self._pharmacophores)

    def __getitem__(self, frame):
        """
            Returns
            -------
            Pharmacophore
        """
        return self._pharmacophores[frame]
