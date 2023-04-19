from collections import namedtuple
import os
from rdkit import RDConfig, RDLogger
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures, rdDistGeom
from rdkit.Chem.Pharm3D import EmbedLib

from openpharmacophore import LigandBasedPharmacophore, \
    LigandReceptorPharmacophore, Pharmacophore, Ligand
from openpharmacophore.pharmacophore.rdkit_pharmacophore import pharmacophore_to_rdkit
from openpharmacophore.pharmacophore.align import align_ligand_to_pharmacophore


Match = namedtuple("Match", ["ligand", "rmsd"])
RDLogger.DisableLog('rdApp.*')  # Disable rdkit warnings


class VirtualScreening:
    """ Class to perform virtual screening with pharmacophores."""

    feature_factory = ChemicalFeatures.BuildFeatureFactory(
        os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))

    def __init__(self, pharmacophores):
        self._pharmacophores = self.add_pharmacophores(pharmacophores)
        self._matches = [[] for _ in range(len(pharmacophores))]
        self._fails = [0] * len(pharmacophores)

    @property
    def n_pharmacophores(self) -> int:
        return len(self._pharmacophores)

    @property
    def pharmacophores(self):
        return self._pharmacophores

    @property
    def fails(self):
        return self._fails

    @property
    def matches(self):
        return self._matches

    @staticmethod
    def add_pharmacophores(pharmacophores):
        """ Add new pharmacophores. """
        pharma_list = []
        for pharma in pharmacophores:
            if isinstance(pharma, (LigandReceptorPharmacophore, LigandBasedPharmacophore)):
                for ph in pharma:
                    pharma_list.append(ph)
            elif isinstance(pharma, Pharmacophore):
                pharma_list.append(pharma)
            else:
                raise TypeError(f"Unexpected type {type(pharma)}")
        return pharma_list

    def screen(self, db):
        """ Screen a database of molecules with each of the pharmacophores.

            Parameters
            ---------
            db : InMemoryMolDB or MolDB
                The database with the molecules.
        """
        # TODO: add progress bar
        for ii, pharma in enumerate(self._pharmacophores):
            pharma = pharmacophore_to_rdkit(pharma)
            for lig in db:
                lig = lig.to_rdkit()
                mappings = self._feature_mappings(lig, pharma)
                if len(mappings) > 0:
                    atom_match = [list(x.GetAtomIds()) for x in mappings]
                    lig = Chem.AddHs(lig)
                    result = align_ligand_to_pharmacophore(lig, atom_match, pharma)
                    if result is not None:
                        aligned_mol, rmsd = result
                        self._matches[ii].append(Match(Ligand(aligned_mol), rmsd))
                    else:
                        self._fails[ii] += 1
                else:
                    self._fails[ii] += 1

    @staticmethod
    def _feature_mappings(ligand, pharmacophore):
        """ Maps the chemical features in the given ligand to those of
            the pharmacophore.

            Parameters
            ----------
            ligand : rdkit.Mol
            pharmacophore : rdkit.Chem.Pharm3D.Pharmacophore

            Returns
            -------
            list [rdkit.MolChemicalFeature]

        """
        can_match, all_matches = EmbedLib.MatchPharmacophoreToMol(
            ligand, VirtualScreening.feature_factory, pharmacophore)
        if can_match:
            bounds_matrix = rdDistGeom.GetMoleculeBoundsMatrix(ligand)
            failed, _, matched_features, _ = EmbedLib.MatchPharmacophore(
                all_matches, bounds_matrix, pharmacophore, useDownsampling=True)
            return matched_features
        return []
