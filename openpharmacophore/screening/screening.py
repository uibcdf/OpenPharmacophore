from openpharmacophore import LigandBasedPharmacophore, \
    LigandReceptorPharmacophore, Pharmacophore


class VirtualScreening:
    """ Class to perform virtual screening with pharmacophores."""

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
        for ii, pharma in enumerate(self._pharmacophores):
            fail = False
            for lig in db:
                lig_feats = lig.get_chem_feats(0)
                if pharma.can_match(lig):
                    mappings = pharma.feature_mappings(lig_feats)
                    if len(mappings) > 0:
                        lig.add_hydrogens()
                        rmsd = pharma.align_ligand(lig, lig_feats)
                        if rmsd is not None:
                            self._matches[ii].append((lig, rmsd))
                        else:
                            fail = True
                    else:
                        fail = True
                else:
                    fail = True

                if fail:
                    self._fails[ii] += 1
