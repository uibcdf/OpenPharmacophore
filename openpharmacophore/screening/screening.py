from ..pharmacophore import LigandBasedPharmacophore, StructureBasedPharmacophore
from rdkit import RDConfig, Geometry
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures, rdDistGeom, rdMolTransforms
from rdkit.Chem.Pharm3D import EmbedLib
from rdkit.Numerics import rdAlignment
from operator import itemgetter
import os


class VirtualScreening:
    """ Class for virtual screening with a pharmacophore.
    """

    feature_factory = ChemicalFeatures.BuildFeatureFactory(
        os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))

    def __init__(self, pharmacophore):
        self._matches = []  # Nested list of molecules that match a pharmacophore
        self._pharmacophores = []  # List of rdkit pharmacophores

        if isinstance(pharmacophore, LigandBasedPharmacophore):
            self._matches.append([])
            self._pharmacophores.append(pharmacophore.to_rdkit())
        elif isinstance(pharmacophore, StructureBasedPharmacophore):
            for ii in range(pharmacophore.num_frames):
                self._matches.append([])
                self._pharmacophores.append(pharmacophore.to_rdkit(ii))

    @property
    def pharmacophores(self):
        return self._pharmacophores

    @property
    def matches(self):
        return self._matches

    @staticmethod
    def _align_to_pharmacophore(molecule, pharmacophore):
        """ Align a molecule to a given pharmacophore.

            Parameters
            ----------
            molecule : rdkit.Mol
                Molecule to align.

            pharmacophore : rdkit.Chem.Pharm3D.Pharmacophore
                An rdkit pharmacophore
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
        _, embeddings, _ = EmbedLib.EmbedPharmacophore(mol_H, atom_match, pharmacophore, count=10)
        SSDs = VirtualScreening._transform_embeddings(pharmacophore, embeddings, atom_match)
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
        align_probe = []   # probe points to align to reference points
        for match_ids in atom_match:
            dummy_point = Geometry.Point3D(0.0, 0.0, 0.0)
            for id_ in match_ids:
                dummy_point += conformer_embed.GetAtomPosition(id_)
            dummy_point /= len(match_ids)
            align_probe.append(dummy_point)

        ssd, transform_matrix = rdAlignment.GetAlignmentTransform(align_ref, align_probe)
        return ssd, transform_matrix
