from openpharmacophore.screening.exceptions import NotAPharmacophoreError
from openpharmacophore import LigandBasedPharmacophore, LigandReceptorPharmacophore
from openpharmacophore.io import mol_file_iterator
from rdkit import RDConfig, Geometry, RDLogger
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures, rdDistGeom, rdMolTransforms
from rdkit.Chem.Pharm3D import EmbedLib
from rdkit.Numerics import rdAlignment
from tqdm.auto import tqdm
from collections import namedtuple
from operator import itemgetter
import os


Match = namedtuple("Match", ["mol", "score"])
RDLogger.DisableLog('rdApp.*')  # Disable rdkit warnings


class VirtualScreening:
    """ Class for virtual screening with a pharmacophore.
    """

    feature_factory = ChemicalFeatures.BuildFeatureFactory(
        os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))

    def __init__(self, pharmacophore):
        self._matches = []  # Nested list of molecules that match a pharmacophore
        self._pharmacophores = []  # List of rdkit pharmacophores
        self._fails = []

        for ii in range(len(pharmacophore)):
            self._matches.append([])
            self._pharmacophores.append(pharmacophore.to_rdkit(ii))
            self._fails.append(0)

        if not isinstance(pharmacophore,
                          (LigandReceptorPharmacophore,
                           LigandBasedPharmacophore)):
            raise NotAPharmacophoreError(type(pharmacophore))

    @property
    def pharmacophores(self):
        return self._pharmacophores

    @property
    def matches(self):
        return self._matches

    def fails(self, index):
        """ Returns the number of molecules that failed to match the pharmacophore
            with given index.

            Parameters
            ----------
            index : int
                The index of the pharmacophore
        """
        return self._fails[index]

    def num_mols(self, index):
        """ Returns the number of molecules that were screened with the pharmacophore
            with given index.

            Parameters
            ----------
            index : int
                The index of the pharmacophore
        """
        return self._fails[index] + len(self._matches[index])

    def from_list(self, molecules, pharmacophore_index):
        """ Screen molecules from a list.

            Parameters
            ----------
            molecules : list[rdkit.Mol]
                The list of molecules that will be screened.

            pharmacophore_index : int
                The index of the pharmacophore that will be used.
        """
        self._from_iterable(molecules, pharmacophore_index)

    def from_file(self, file_name, pharmacophore_index):
        """ Screen molecules from a file.

            Parameters
            ----------
            file_name : str
                Name of the file

            pharmacophore_index : int
                The index of the pharmacophore that will be used.
        """
        molecules = mol_file_iterator(file_name)
        self._from_iterable(molecules, pharmacophore_index)

    def from_dir(self, path, pharmacophore_index, skip=None):
        """ Screen molecules from a directory.

            Parameters
            ----------
            path : str
                Name of the directory.

            pharmacophore_index : int
                The index of the pharmacophore that will be used.

            skip : list[str]
                Skip files in this list.
        """
        # TODO: add progress bar
        file_formats = ["smi", "mol2", "sdf"]
        for root, dirs, files in os.walk(path):
            for file in files:
                if file in skip:
                    continue
                if file.split(".")[-1] in file_formats:
                    self.from_file(os.path.join(root, file), pharmacophore_index)

    def _from_iterable(self, molecules_iter, pharmacophore_index):
        """ Screen molecules from an iterator.

           Parameters
           ----------
           molecules_iter : Iterator[rdkit.Mol]
               An iterator of molecules.

           pharmacophore_index : int
               The index of the pharmacophore that will be used.
       """

        for mol in tqdm(molecules_iter):
            if mol is not None:
                mol_and_score = self._align_to_pharmacophore(
                    mol, self._pharmacophores[pharmacophore_index])
                if mol_and_score is not None:
                    match = Match(mol_and_score[0], mol_and_score[1])
                    self._matches[pharmacophore_index].append(match)
                else:
                    self._fails[pharmacophore_index] += 1

    @staticmethod
    def _align_to_pharmacophore(molecule, pharmacophore):
        """ Align a molecule to a given pharmacophore. If the molecule
            can be aligned it returns the molecule with a conformer that
            matches the pharmacophore and the SSD value of the alignment.

            If the molecule can't be matched to the pharmacophore, None is
            returned.

            Parameters
            ----------
            molecule : rdkit.Mol
                Molecule to align.

            pharmacophore : rdkit.Chem.Pharm3D.Pharmacophore
                An rdkit pharmacophore.

            Returns
            -------
            rdkit.Mol
                The molecule with a conformer that matches the pharmacophore.

            float
                The SSD value of the alignment.
        """
        bounds_matrix = rdDistGeom.GetMoleculeBoundsMatrix(molecule)
        # Check if the molecule features can match with the pharmacophore.
        can_match, all_matches = EmbedLib.MatchPharmacophoreToMol(
            molecule, VirtualScreening.feature_factory, pharmacophore)
        # all_matches is a list of tuples where each tuple contains the chemical features
        if can_match:
            # Match the molecule to the pharmacophore without aligning it
            failed, _, matched_features, _ = EmbedLib.MatchPharmacophore(
                all_matches, bounds_matrix, pharmacophore, useDownsampling=True)
            if failed:
                return
        else:
            return
        atom_match = [list(x.GetAtomIds()) for x in matched_features]
        mol_H = Chem.AddHs(molecule)
        # Embed molecule onto the pharmacophore
        # embeddings is a list of molecules with a single conformer
        try:
            _, embeddings, _ = EmbedLib.EmbedPharmacophore(mol_H, atom_match, pharmacophore, count=10)
        except KeyError:
            # When embed fails it raises a key error
            return
        SSDs = VirtualScreening._transform_embeddings(pharmacophore, embeddings, atom_match)
        if len(SSDs) == 0:
            return
        best_fit_index = min(enumerate(SSDs), key=itemgetter(1))[0]

        return embeddings[best_fit_index], SSDs[best_fit_index]

    @staticmethod
    def _transform_embeddings(pharmacophore, embeddings, atom_match):
        """ Transform embeddings. Performs the alignment of the conformers of
            a molecule to the pharmacophore.

            Parameters
            ----------
            pharmacophore: rdkit.Chem.Pharm3D.Pharmacophore
                A pharmacophore object.

            embeddings: list[rdkit.Mol]
                List of molecules (the same molecule) with a single conformer.

            atom_match: list[list[int]]
                A nested list of atoms ids that match the pharmacophore.

            Returns
            -------
            SSDs: list[float]
                List of sum of square deviations (SSD) values for the alignments.

        """
        align_ref = [f.GetPos() for f in pharmacophore.getFeatures()]
        ssds = []
        for embedding in embeddings:
            conformer = embedding.GetConformer()
            ssd, transform_matrix = VirtualScreening._transform_matrix(align_ref, conformer, atom_match)
            # Transform the coordinates of the conformer
            rdMolTransforms.TransformConformer(conformer, transform_matrix)
            ssds.append(ssd)

        return ssds

    @staticmethod
    def _transform_matrix(align_ref, conformer_embed, atom_match):
        """ Get the transformation matrix of a conformer that is aligned to
            a pharmacophore.

            Parameters
            ----------
            align_ref: list[rdkit.Geometry.Point3D]
                list of pharmacophore reference points for the alignment.

            conformer_embed: rdkit.Conformer
                The conformer embedding.

            atom_match: list[list]
                Nested list of atoms ids that match the pharmacophore.

            Returns
            -------
            ssd: float
                SSD (sum of square deviations) value of the alignment.

            transform_matrix: numpy.ndarray; shape(4, 4)
                The transform matrix.
        """
        align_probe = []  # probe points to align to reference points
        for match_ids in atom_match:
            dummy_point = Geometry.Point3D(0.0, 0.0, 0.0)
            for id_ in match_ids:
                dummy_point += conformer_embed.GetAtomPosition(id_)
            dummy_point /= len(match_ids)
            align_probe.append(dummy_point)

        ssd, transform_matrix = rdAlignment.GetAlignmentTransform(align_ref, align_probe)
        return ssd, transform_matrix
