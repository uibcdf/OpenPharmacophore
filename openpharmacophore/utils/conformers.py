import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "3-clause BSD"

# options
FORCEFIELD_CHOICES = ("uff", "mmff94", "mmff94s")


class ConformerGenerator(object):
    """Generate conformers using RDKit.

    Procedure
    ---------
    1. Generate a pool of conformers.
    2. Minimize conformers.
    3. Filter conformers using an RMSD threshold and optional minimum energy
       difference.
    Note that pruning is done _after_ minimization, which differs from the
    protocol described in the references.
    References
    ----------
    * https://github.com/keiserlab/e3fp/blob/master/e3fp/conformer/generator.py
    * http://rdkit.org/docs/GettingStartedInPython.html
      #working-with-3d-molecules
    * http://pubs.acs.org/doi/full/10.1021/ci2004658
    * https://github.com/skearnes/rdkit-utils/blob/master/rdkit_utils/
      conformers.py
    """

    def __init__(
        self,
        num_conf=-1,
        rmsd_cutoff=0.5,
        max_energy_diff=-1.0,
        forcefield="uff",
        pool_multiplier=1,
        seed=-1,
    ):
        """Initialize generator settings.

        Parameters
        ----------
        num_conf : int, optional
            Maximum number of conformers to generate (after pruning). -1
            results in auto selection of max_conformers.
        pool_multiplier : int, optional
            Factor to multiply by max_conformers to generate the initial
            conformer pool. Since conformers are filtered after energy
            minimization, increasing the size of the pool increases the chance
            of identifying max_conformers unique conformers.
        rmsd_cutoff : float, optional
            RMSD cutoff for pruning conformers. If None or negative, no
            pruning is performed.
        max_energy_diff : float, optional
            If set, conformers with energies this amount above the minimum
            energy conformer are not accepted.
        forcefield : {'uff', 'mmff94', 'mmff94s'}, optional
            Force field to use for conformer energy calculation and
            minimization.
        seed : int, optional
            Random seed for conformer generation. If -1, the random number
            generator is unseeded.
        """
        self.max_conformers = num_conf
        if not rmsd_cutoff or rmsd_cutoff < 0:
            rmsd_cutoff = -1.0
        self.rmsd_cutoff = rmsd_cutoff

        if max_energy_diff is None or max_energy_diff < 0:
            max_energy_diff = -1.0
        self.max_energy_diff = max_energy_diff

        if forcefield not in FORCEFIELD_CHOICES:
            raise ValueError(
                "%s is not a valid option for forcefield" % forcefield
            )
        self.forcefield = forcefield
        self.pool_multiplier = pool_multiplier
        self.seed = seed

    def __call__(self, mol):
        """Generate conformers for a molecule.
        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        Returns
        -------
        RDKit Mol : copy of the input molecule with embedded conformers
        """
        return self.generate_conformers(mol)

    def generate_conformers(self, mol):
        """ Generate conformers for a molecule.

        Parameters
        ----------
        mol : rdkit.Chem.Mol
            Molecule.
        Returns
        -------
        rdkit.Chem.Mol
            Molecule with embedded conformers
        """
        # initial embedding
        mol = self.embed_molecule(mol)
        if mol.GetNumConformers() == 0:
            raise ValueError("Failed to generate conformers")

        # minimization and filtering
        self.minimize_conformers(mol)
        mol, indices, energies, rmsds = self.filter_conformers(mol)
        mol = Chem.AddHs(mol, addCoords=True)
        return mol

    @staticmethod
    def get_num_conformers(mol):
        """Return ideal number of conformers from rotatable bond number in model.

        Parameters
        ----------
        mol : rdkit.Chem.Mol
            RDKit `Mol` object for molecule

        Returns
        ------
        num_conf : int
            Target number of conformers to accept
        """
        num_rot = AllChem.CalcNumRotatableBonds(mol)
        if num_rot < 8:
            return 50
        elif 8 <= num_rot <= 12:
            return 200
        elif num_rot > 12:
            return 300
        else:
            return 0

    def embed_molecule(self, mol):
        """ Generate conformers, possibly with pruning.

            Parameters
            ----------
            mol : RDKit Mol
                Molecule.
        """
        if not self.has_hydrogens(mol):
            mol = Chem.AddHs(mol)  # add hydrogens
        Chem.SanitizeMol(mol)
        if self.max_conformers == -1 or type(self.max_conformers) is not int:
            self.max_conformers = self.get_num_conformers(mol)
        n_confs = self.max_conformers * self.pool_multiplier
        AllChem.EmbedMultipleConfs(
            mol,
            numConfs=n_confs,
            maxAttempts=10 * n_confs,
            pruneRmsThresh=-1.0,
            randomSeed=self.seed,
            ignoreSmoothingFailures=True,
        )
        return mol

    def get_molecule_force_field(self, mol, conf_id=None, **kwargs):
        """Get a force field for a molecule.
        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        conf_id : int, optional
            ID of the conformer to associate with the force field.
        **kwargs : dict, optional
            Keyword arguments for force field constructor.
        """
        if self.forcefield == "uff":
            ff = AllChem.UFFGetMoleculeForceField(
                mol, confId=conf_id, **kwargs
            )
        elif self.forcefield.startswith("mmff"):
            AllChem.MMFFSanitizeMolecule(mol)
            mmff_props = AllChem.MMFFGetMoleculeProperties(
                mol, mmffVariant=self.forcefield
            )
            ff = AllChem.MMFFGetMoleculeForceField(
                mol, mmff_props, confId=conf_id, **kwargs
            )
        else:
            raise ValueError(
                "Invalid forcefield " + "'{}'.".format(self.forcefield)
            )
        return ff

    def minimize_conformers(self, mol):
        """ Minimize molecule conformers.
        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        for conf in mol.GetConformers():
            ff = self.get_molecule_force_field(mol, conf_id=conf.GetId())
            ff.Minimize()

    def get_conformer_energies(self, mol):
        """ Calculate conformer energies.
        Parameters
        ----------
        mol : rdkit.Chem.Mol
            Molecule.

        Returns
        -------
        energies : array_like
            Minimized conformer energies.
        """
        num_conf = mol.GetNumConformers()
        energies = np.empty((num_conf,), dtype=float)
        for i, conf in enumerate(mol.GetConformers()):
            ff = self.get_molecule_force_field(mol, conf_id=conf.GetId())
            energies[i] = ff.CalcEnergy()
        return energies

    def filter_conformers(self, mol):
        """ Filter conformers which do not meet an RMSD threshold.

            Parameters
            ----------
            mol : rdkit.Chem.Mol
                Molecule.
            Returns
            -------
            rdkit.Chem.Mol
                A new RDKit Mol containing the chosen conformers, sorted by
                increasing energy.
        """
        energies = self.get_conformer_energies(mol)
        energy_below_threshold = np.ones_like(energies, dtype=np.bool_)

        sort = np.argsort(energies)  # sort by increasing energy
        confs = np.array(mol.GetConformers())

        mol = Chem.RemoveHs(mol)
        accepted = []  # always accept lowest-energy conformer
        rejected = []
        rmsds = np.zeros((confs.shape[0], confs.shape[0]), dtype=float)
        for i, fit_ind in enumerate(sort):
            accepted_num = len(accepted)

            # always accept lowest-energy conformer
            if accepted_num == 0:
                accepted.append(fit_ind)

                # pre-compute if Es are in acceptable range of min E
                if self.max_energy_diff != -1.0:
                    energy_below_threshold = (
                        energies <= energies[fit_ind] + self.max_energy_diff
                    )

                continue

            # check if energy is too high
            if not energy_below_threshold[fit_ind]:
                rejected.append(fit_ind)
                continue

            # get RMSD to selected conformers
            these_rmsds = np.zeros((accepted_num,), dtype=float)
            # reverse so all confs aligned to lowest energy
            for j, accepted_ind in self.reverse_enumerate(accepted):
                this_rmsd = AllChem.GetBestRMS(
                    mol,
                    mol,
                    confs[accepted_ind].GetId(),
                    confs[fit_ind].GetId(),
                )
                # reject conformers within the RMSD threshold
                if this_rmsd < self.rmsd_cutoff:
                    rejected.append(fit_ind)
                    break
                else:
                    these_rmsds[-j - 1] = this_rmsd
            else:
                rmsds[fit_ind, accepted] = these_rmsds
                rmsds[accepted, fit_ind] = these_rmsds
                accepted.append(fit_ind)

        # slice and order rmsds and energies to match accepted list
        rmsds = rmsds[np.ix_(accepted, accepted)]
        energies = energies[accepted]

        # create a new molecule with all conformers, sorted by energy
        new = Chem.Mol(mol)
        new.RemoveAllConformers()
        conf_ids = [conf.GetId() for conf in mol.GetConformers()]
        for i in accepted:
            conf = mol.GetConformer(conf_ids[i])
            new.AddConformer(conf, assignId=True)

        return new, np.asarray(accepted, dtype=int), energies, rmsds

    @staticmethod
    def reverse_enumerate(iterable):
        """Enumerate, but with the last result first but still numbered last.
        Parameters
        ----------
        iterable : some 1-D iterable
        Returns
        -------
        iterable:
            Reverse of `enumerate` function
        """
        return zip(reversed(range(len(iterable))), reversed(iterable))

    @staticmethod
    def has_hydrogens(mol):
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "H":
                return True
            return False
