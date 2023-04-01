# OpenPharmacophore

[![GitHub Actions Build Status](https://github.com/uibcdf/OpenPharmacophore/workflows/CI/badge.svg)](https://github.com/uibcdf/molsysmt/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/uibcdf/OpenPharmacophore/master/graph/badge.svg)](https://codecov.io/gh/uibcdf/OpenPharmacophore/branch/master)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


OpenPharmacophore is a library in progress to work with pharmacophore models and virtual screening. 
Its main goal is to derive pharmacophores from ligand, ligand-receptor, and receptor based methods.
It will also offer the possibility to derive pharmacophores from molecular dynamics trajectories.


<img src="./pharmacophore.png" alt="pharmacophore">


## Installation


The latest "stable" version of OpenPharmacophore can be installed from the UIBCDF Anaconda channel:

```bash
conda -c uibcdf openpharmacophore
```

If you want to work with the not so tested last beta version, the installation command is the following:

```bash
conda install -c uibcdf/label/dev openpharmacophore
```

The former beta version is nothing but a quenched version from the main github repository of this project which it is done from time to time with few scruples. The raw code fully alive can be installed from this github repo as follows:

```bash
git clone https://github.com/uibcdf/openpharmacophore.git
cd OpenPharmacophore
python setup.py develop
```

In the first two cases, OpenPharmacophore can be uninstalled with conda:

```bash
conda remove openpharmacophore
```

But if you installed it straight from its github central repository, do the following to uninstall it:

```bash
pip uninstall openpharmacophore
```


## License


## Acknowledgments and Copyright

### Copyright

Copyright (c) 2021 to 2023, UIBCDF

## Contributors

- Daniel Ibarrola Sánchez (UAM-I/HIMFG)    
- Rafael A. Zubillaga Luna (UAM-I)    
- Liliana M. Moreno Vargas (HIMFG)    
- Diego Prada Gracia (HIMFG)    

#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
