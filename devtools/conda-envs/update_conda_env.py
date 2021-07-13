import argparse
import os
import re
import glob
import shutil
import subprocess as sp
from tempfile import TemporaryDirectory
from contextlib import contextmanager


# Args
parser = argparse.ArgumentParser(description='Updates the activated conda environment with the packages in a yaml file')
parser.add_argument('conda_file',
                    help='The file for the created Python environment')

args = parser.parse_args()

# Figure out conda path
if "CONDA_EXE" in os.environ:
    conda_path = os.environ["CONDA_EXE"]
else:
    conda_path = shutil.which("conda")
if conda_path is None:
    raise RuntimeError("Could not find a conda binary in CONDA_EXE variable or in executable search path")

print("CONDA FILE NAME {}".format(args.conda_file))
print("CONDA PATH      {}".format(conda_path))

# Write to a temp directory which will always be cleaned up

sp.call("{} env update --file {} --prune".format(conda_path, args.conda_file), shell=True)

