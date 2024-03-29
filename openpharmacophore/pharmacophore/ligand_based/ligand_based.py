from openpharmacophore import Ligand, Pharmacophore
import openpharmacophore.pharmacophore.ligand_based.common_pharmacophore as cp


class LigandBasedPharmacophore:
    """ Class to store and extract pharmacophores from a set of ligands.

        Parameters
        ----------
        ligands : list[Ligand]

    """
    # TODO: add exclusion volumes to ligand based model
    def __init__(self, ligands):
        self._ligands = ligands
        self._pharmacophores = []  # type: list[Pharmacophore]

    @property
    def pharmacophores(self):
        return self._pharmacophores

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

        # TODO: use chemical features directionality
        chem_feats = self._get_chem_feats(self._ligands)
        extractor = cp.CommonPharmacophoreFinder(*args, **kwargs)
        self._pharmacophores = extractor(chem_feats, n_points, min_actives, max_pharmacophores)

    @staticmethod
    def _get_chem_feats(ligands):
        """ Get the chemical feature of all ligands in the pharmacophore.

            Parameters
            ----------
            list[Ligand]

            Returns
            -------
            list[ChemFeatContainer]
        """
        chem_feats = []
        for ii, lig in enumerate(ligands):
            ligand_feats = []
            for conf in range(lig.n_conformers):
                chem_f = lig.get_chem_feats(conf, indices=False)
                ligand_feats.append(chem_f)
            chem_feats.append(ligand_feats)
        return chem_feats

    def __len__(self):
        return len(self._pharmacophores)

    def __getitem__(self, frame):
        """
            Returns
            -------
            Pharmacophore
        """
        return self._pharmacophores[frame]
