# Installation

The latest "stable" version of Molecular Systems can be installed from the UIBCDF Anaconda channel:

```bash
conda -c uibcdf molecular_systems
```

If you want to work with the not so tested last beta version, the installation command is the following:

```bash
conda install -c uibcdf/label/dev molecular_systems
```

The former beta version is nothing but a quenched version from the main github repository of this project which it is done from time to time with few scruples. The raw code fully alive can be installed from this github repo as follows:

```bash
git clone https://github.com/uibcdf/Molecular_Systems.git
cd Molecular_Systems
python setup.py develop
```

In the first two cases, Molecular Systems can be uninstalled with conda:

```bash
conda remove molecular_systems
```

But if you installed Molecular Systems straight from its github central repository, do the following to uninstall it:

```bash
pip uninstall molecular_systems
```
