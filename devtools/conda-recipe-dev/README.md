# Instructions

## Required

```bash
conda install conda-build anaconda-client
```

## Building and pushing to https://anaconda.org/uibcdf

```bash
conda build . --no-anaconda-upload
PACKAGE_OUTPUT=`conda build . --output`
anaconda login
anaconda upload --user uibcdf $PACKAGE_OUTPUT --label dev
conda build purge
anaconda logout
```

## Install

```
conda install -c uibcdf/dev openpharmacophore
```

## Additional Info
https://docs.anaconda.com/anaconda-cloud/user-guide/tasks/work-with-packages
