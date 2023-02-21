from mdtraj.formats.pdb.pdbfile import PDBTrajectoryFile, _format_83
import pyunitwizard as puw
import numpy as np
import io
from typing import Union


class PDBFileWriter(PDBTrajectoryFile):
    """ Class to write pdb files.

        Inherits from mdtraj.PDBTrajectory file. Fixes bug with atom positions.
        Also allows writing to a stringio.

    """

    def __init__(self, file: Union[io.FileIO, io.StringIO]):
        self._open = False
        self._file = None
        self._topology = None
        self._positions = None
        self._mode = "w"
        self._last_topology = None
        self._standard_names = True

        self._header_written = False
        self._footer_written = False

        self._file = file

        self._open = True

    def write(self, positions, topology, modelIndex=None, unitcell_lengths=None,
              unitcell_angles=None, bfactors=None):
        """Write a PDB file to disk
        Parameters
        ----------
        positions : array_like
            The list of atomic positions to write. Shape (n_atoms, 3)
        topology : mdtraj.Topology
            The Topology defining the model to write.
        modelIndex : {int, None}
            If not None, the model will be surrounded by MODEL/ENDMDL records
            with this index
        unitcell_lengths : {tuple, None}
            Lengths of the three unit cell vectors, or None for a non-periodic system
        unitcell_angles : {tuple, None}
            Angles between the three unit cell vectors, or None for a non-periodic system
        bfactors : array_like, default=None, shape=(n_atoms,)
            Save bfactors with pdb file. Should contain a single number for
            each atom in the topology
        """
        if not self._mode == 'w':
            raise ValueError('file not opened for writing')
        if not self._header_written:
            self._write_header(unitcell_lengths, unitcell_angles)
            self._header_written = True

        if topology.n_atoms != positions.shape[0]:
            raise ValueError(
                f'The number of positions must match the number of atoms.'
                f'Num atoms: {topology.n_atoms}. Positions shape: {positions.shape}'
            )
        if np.any(np.isnan(positions)):
            raise ValueError('Particle position is NaN')
        if np.any(np.isinf(positions)):
            raise ValueError('Particle position is infinite')

        self._last_topology = topology

        if bfactors is None:
            bfactors = ['{0:5.2f}'.format(0.0)] * len(positions)
        else:
            if (np.max(bfactors) >= 100) or (np.min(bfactors) <= -10):
                raise ValueError("bfactors must be in (-10, 100)")

            bfactors = ['{0:5.2f}'.format(b) for b in bfactors]

        atomIndex = 1
        if modelIndex is not None:
            print("MODEL     %4d" % modelIndex, file=self._file)
        for (chainIndex, chain) in enumerate(topology.chains):
            chainName = self._chain_names[chainIndex % len(self._chain_names)]
            residues = list(chain.residues)
            for (resIndex, res) in enumerate(residues):
                if len(res.name) > 3:
                    resName = res.name[:3]
                else:
                    resName = res.name
                for atom in res.atoms:
                    if len(atom.name) < 4 and atom.name[:1].isalpha() and (
                            atom.element is None or len(atom.element.symbol) < 2):
                        atomName = ' ' + atom.name
                    elif len(atom.name) > 4:
                        atomName = atom.name[:4]
                    else:
                        atomName = atom.name
                    coords = positions[atom.index]
                    if atom.element is not None:
                        symbol = atom.element.symbol
                    else:
                        symbol = ' '
                    if atom.serial is not None and len(topology._chains) < 2:
                        atomSerial = atom.serial
                    else:
                        atomSerial = atomIndex
                    line = "ATOM  %5d %-4s %3s %1s%4d    %s%s%s  1.00 %5s      %-4s%2s  " % (
                        # Right-justify atom symbol
                        atomSerial % 100000, atomName, resName, chainName,
                        (res.resSeq % 10000), _format_83(coords[0]),
                        _format_83(coords[1]), _format_83(coords[2]),
                        bfactors[atom.index], atom.segment_id[:4], symbol[-2:]
                    )
                    assert len(line) == 80, 'Fixed width overflow detected'
                    print(line, file=self._file)
                    atomIndex += 1
                if resIndex == len(residues) - 1:
                    print("TER   %5d      %3s %s%4d" % (atomSerial + 1, resName, chainName, res.resSeq),
                          file=self._file)
                    atomIndex += 1

        if modelIndex is not None:
            print("ENDMDL", file=self._file)

    @property
    def file(self) -> Union[io.FileIO, io.StringIO]:
        return self._file


def write_pdb_block(topology, coords):
    """ Write a pdb block for a molecule.

        Parameters
        ----------
        topology : Topology
            The topology of the molecule.

        coords : QuantityLike
            The positions of the atoms.
    """
    with PDBFileWriter(file=io.StringIO()) as pdb:
        pdb.write(puw.get_value(coords[0], "angstroms"),
                  topology.top,
                  modelIndex=0,
                  bfactors=None
                  )

        pdb_block = pdb.file.getvalue()

    return pdb_block


def write_pdb(file_name, topology, coords):
    """ Write a pdb block for a molecule.

            Parameters
            ----------
            file_name: str
                Name of the file to write.

            topology : Topology
                The topology of the molecule.

            coords : QuantityLike
                The positions of the atoms.
        """
    with PDBFileWriter(file=open(file_name, "w")) as pdb:
        pdb.write(puw.get_value(coords[0], "angstroms"),
                  topology.top,
                  modelIndex=0,
                  bfactors=None
                  )
