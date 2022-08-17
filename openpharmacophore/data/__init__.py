# Loads all data files in this folder and puts them into
# appropriate dictionaries.

from pathlib import PurePath
import os

parent = PurePath(__file__).parent

bioassays = {}
ligands = {}
pharmacophores = {}
pdb = {}
trajectories = {}
topologies = {}
smarts_features = parent.joinpath("smarts_features.txt")

# Traverse all folders from the data directory and
# fill the dictionaries
for dir_path, _, filenames in os.walk(parent):

    for file_name in filenames:
        name = file_name.split(".")[0]
        if "bioassays" in dir_path:
            bioassays[name] = str(parent.joinpath(
                os.path.join(dir_path, file_name)
            ))
        elif "pharmacophores" in dir_path:
            pharmacophores[name] = str(parent.joinpath(
                os.path.join(dir_path, file_name)
            ))
        elif "ligands" in dir_path:
            ligands[name] = str(parent.joinpath(
                os.path.join(dir_path, file_name)
            ))
        elif "pdb" in dir_path:
            pdb[name] = str(parent.joinpath(
                os.path.join(dir_path, file_name)
            ))
        elif "trajectories" in dir_path:
            trajectories[name] = str(parent.joinpath(
                os.path.join(dir_path, file_name)
            ))
        elif "topologies" in dir_path:
            topologies[name] = str(parent.joinpath(
                os.path.join(dir_path, file_name)
            ))
