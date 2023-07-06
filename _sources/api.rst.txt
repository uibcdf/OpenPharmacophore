API Documentation
=================

Molecular Systems
-----------------------

.. autosummary::
   :toctree: _autosummary
   :template: custom-class-template.rst
   :recursive:

   openpharmacophore.Protein
   openpharmacophore.Ligand
   openpharmacophore.ComplexBindingSite

Pharmacophores
---------------
.. autosummary::
   :toctree: _autosummary
   :template: custom-class-template.rst
   :recursive:

    openpharmacophore.PharmacophoricPoint
    openpharmacophore.Pharmacophore
    openpharmacophore.LigandBasedPharmacophore
    openpharmacophore.LigandReceptorPharmacophore

Visualization
--------------
.. autosummary::
   :toctree: _autosummary
   :template: custom-class-template.rst
   :recursive:

    openpharmacophore.Viewer

.. autosummary::
   :toctree: _autosummary
   :template: base.rst
   :recursive:

    openpharmacophore.draw_ligands
    openpharmacophore.draw_ligands_chem_feats

File IO
-------
.. autosummary::
   :toctree: _autosummary
   :template: base.rst
   :recursive:

    openpharmacophore.io.pharmacophore_reader.read_mol2
    openpharmacophore.io.pharmacophore_reader.read_ph4
    openpharmacophore.io.pharmacophore_reader.read_json
    openpharmacophore.io.pharmacophore_reader.read_ligandscout

    openpharmacophore.io.pharmacophore_writer.save_mol2
    openpharmacophore.io.pharmacophore_writer.save_ph4
    openpharmacophore.io.pharmacophore_writer.save_json
    openpharmacophore.io.pharmacophore_writer.save_pml

    openpharmacophore.io.mol_files.read_mol2
    openpharmacophore.io.mol_files.mol2_supplier
    openpharmacophore.io.mol_files.read_sdf


Virtual Screening
-----------------
.. autosummary::
   :toctree: _autosummary
   :template: custom-class-template.rst
   :recursive:

    openpharmacophore.VirtualScreening
    openpharmacophore.screening.mol_db.InMemoryMolDB
    openpharmacophore.screening.mol_db.MolDB


Others
--------
.. autosummary::
   :toctree: _autosummary
   :template: base.rst
   :recursive:

    openpharmacophore.load
    openpharmacophore.load_ligands
    openpharmacophore.smiles_from_pdb_id
