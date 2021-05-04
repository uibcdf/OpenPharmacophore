from setuptools import find_packages
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension

extensions_list=[]

setup(
    name='openpharmacophore',
    version='0.0.1',
    author='UIBCDF Lab and contributors',
    author_email='uibcdf@gmail.com',
    package_dir={'openpharmacophore': 'openpharmacophore'},
    packages=find_packages(),
    ext_modules=extensions_list,
    package_data={'openpharmacophore': []},
    scripts=[],
    url='http://uibcdf.org',
    download_url ='https://github.com/uibcdf/OpenPharmacophore',
    license='MIT',
    description="---",
    long_description="---",
)

