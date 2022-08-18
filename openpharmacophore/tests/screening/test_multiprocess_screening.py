import openpharmacophore as oph
import openpharmacophore.data as data

# TODO: MultiprocessScreeningClass is not properly tested


def is_substring_in_list(list_, substring):
    for element in list_:
        if substring in element:
            return True
    return False


def test_get_files():
    path_to_files = data.parent.joinpath("ligands")
    file_queue = oph.MultiProcessVirtualScreening._get_files(path_to_files)
    file_list = []
    for ii in range(5):
        file_list.append(file_queue.get())

    assert is_substring_in_list(file_list, "ace.mol2")
    assert is_substring_in_list(file_list, "BAAAML.smi")
    assert is_substring_in_list(file_list, "clique_detection.smi")
    assert is_substring_in_list(file_list, "er_alpha_ligands.sdf")
    assert is_substring_in_list(file_list, "mols.smi")
