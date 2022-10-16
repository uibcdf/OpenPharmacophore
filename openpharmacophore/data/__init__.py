# Loads all data files in this folder and puts them into
# appropriate dictionaries.

from pathlib import PurePath
import os

parent = PurePath(__file__).parent

ligands = {}
pharmacophores = {}
pdb = {}
trajectories = {}
zinc = {}


# Traverse all folders from the data directory and
# fill the dictionaries
for dir_path, _, filenames in os.walk(parent):

    for file_name in filenames:
        if "pharmacophores" in dir_path:
            pharmacophores[file_name] = str(parent.joinpath(
                os.path.join(dir_path, file_name)
            ))
        elif "ligands" in dir_path:
            ligands[file_name] = str(parent.joinpath(
                os.path.join(dir_path, file_name)
            ))
        elif "pdb" in dir_path:
            pdb[file_name] = str(parent.joinpath(
                os.path.join(dir_path, file_name)
            ))
        elif "trajectories" in dir_path:
            trajectories[file_name] = str(parent.joinpath(
                os.path.join(dir_path, file_name)
            ))
        elif "zinc" in dir_path:
            zinc[file_name] = str(parent.joinpath(
                os.path.join(dir_path, file_name)
            ))
