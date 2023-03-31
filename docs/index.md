# OpenPharmacophore


OpenPharmacophore is a library to work with pharmacophore models and perform virtual screening. It 
can derive pharmacophores from ligand-based as well as structured-based methods. It also offers
the possibility to derive pharmacophores from molecular dynamics trajectories.

<img src="_static/pharmacophore.png" alt="pharmacophore">

## Installation


The latest "stable" version of OpenPharmacophore can be installed from the UIBCDF Anaconda channel:

```bash
conda -c uibcdf openpharmacophore
```

## Examples

### Ligand based pharmacophores
  - [Preparing ligand for pharmacophore extraction](contents/examples/ligand-based/ligand_preparation.ipynb)

### Ligand-receptor based pharmacophores
  - [examples of the protein-ligand complex of estrogen receptor with estradiol](contents/examples/ligand-receptor/ligand_receptor_pharmacophore.ipynb)

### Dynamic pharmacophores
  - [Obtaining a ligand-receptor pharmacophore model from a MD trajectory](contents/examples/ligand-receptor/lr_dynamic_pharmacophore.ipynb)

### Receptor based pharmacophores

- Coming soon...

### Virtual Screening

- Coming soon...

### Others

#### Protein-ligand interactions

- [Aromatic Interactions](contents/examples/pl-interactions/aromatic_interactions.ipynb)
- [Hydrogen Bonding](contents/examples/pl-interactions/hydrogen_bonding.ipynb)
- [Hydrophobic Interactions](contents/examples/pl-interactions/hydrophobic_interactions.ipynb)


#### Binding site

- [Protein-ligand complex binding site](contents/examples/binding-site/complex_binding_site.ipynb)


```{eval-rst}
.. toctree::
   :name: examples
   :caption: Examples
   :maxdepth: 2
   :hidden:

   contents/index.md

.. toctree::
   :name: api_doc
   :caption: API Documentation
   :maxdepth: 2
   :hidden:

   api.rst
```
   