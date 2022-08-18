import openpharmacophore.io as io
import openpharmacophore.data as data
import numpy as np
import pyunitwizard as puw
from example_pharmacophores import two_element_pharmacophore, five_element_pharmacophore


def test_from_pharmer():
    points, molecular_system, ligand = io.from_pharmer(data.pharmacophores["1M70"],
                                                       load_mol_sys=False)
    assert len(points) == 5
    assert molecular_system is None
    assert ligand is None

    assert points[0].feature_name == "hb acceptor"
    assert points[0].has_direction
    assert puw.get_value(points[0].radius, "angstroms") == 0.9999999999999999
    assert np.allclose(puw.get_value(points[0].center, "angstroms"),
                       np.array([21.352, -14.531, 19.625]))
    assert np.allclose(points[0].direction,
                       np.array([-0.6405836470264256, 0.7029084735090229, -0.3091476492414897]))

    assert points[1].feature_name == "hb acceptor"
    assert points[1].has_direction
    assert puw.get_value(points[1].radius, "angstroms") == 0.9999999999999999
    assert np.allclose(puw.get_value(points[1].center, "angstroms"),
                       np.array([19.355, -18.32, 23.987]))
    assert np.allclose(points[1].direction,
                       np.array([0.6859059711903811, 0.09092493673854565, 0.721987295292979]))

    assert points[2].feature_name == "hb donor"
    assert points[2].has_direction
    assert puw.get_value(points[2].radius, "angstroms") == 0.9999999999999999
    assert np.allclose(puw.get_value(points[2].center, "angstroms"),
                       np.array([20.977, -16.951, 18.746]))
    assert np.allclose(points[2].direction,
                       np.array([0.71662539652105, -0.5202950182802607, -0.46447942367105033]))

    assert points[3].feature_name == "negative charge"
    assert not points[3].has_direction
    assert puw.get_value(points[3].radius, "angstroms") == 1.5
    assert np.allclose(puw.get_value(points[3].center, "angstroms"),
                       np.array([21.66899, -15.077667, 20.608334]))

    assert points[4].feature_name == "negative charge"
    assert not points[4].has_direction
    assert puw.get_value(points[4].radius, "angstroms") == 1.9999999999999998
    assert np.allclose(puw.get_value(points[4].center, "angstroms"),
                       np.array([19.985, -19.404402, 22.8422]))


def test_to_pharmer_with_two_element_pharmacophore():
    pharmer = io.pharmer._pharmer_dict(two_element_pharmacophore())

    # Expected output from to_pharmer
    expected = {"points": []}

    aromatic_1 = {}
    aromatic_1["name"] = "Aromatic"
    aromatic_1["hasvec"] = True
    aromatic_1["svector"] = {}
    aromatic_1["svector"]["x"] = 0.0
    aromatic_1["svector"]["y"] = 0.0
    aromatic_1["svector"]["z"] = 1.0
    aromatic_1["x"] = 0.9999999999999999
    aromatic_1["y"] = 0.0
    aromatic_1["z"] = 0.0
    aromatic_1["radius"] = 0.9999999999999999
    aromatic_1["enabled"] = True
    aromatic_1["vector_on"] = 0
    aromatic_1["minsize"] = ""
    aromatic_1["maxsize"] = ""
    aromatic_1["selected"] = False

    acceptor = {}
    acceptor["name"] = "HydrogenAcceptor"
    acceptor["hasvec"] = True
    acceptor["svector"] = {}
    acceptor["svector"]["x"] = 0.0
    acceptor["svector"]["y"] = 0.7071067811865475
    acceptor["svector"]["z"] = 0.7071067811865475
    acceptor["x"] = 0.9999999999999999
    acceptor["y"] = 1.9999999999999998
    acceptor["z"] = 1.9999999999999998
    acceptor["radius"] = 0.9999999999999999
    acceptor["enabled"] = True
    acceptor["vector_on"] = 0
    acceptor["minsize"] = ""
    acceptor["maxsize"] = ""
    acceptor["selected"] = False

    expected["points"].append(acceptor)
    expected["points"].append(aromatic_1)
    assert pharmer == expected

    # Test pharmacophore with five elements
    pharmer = io.pharmer._pharmer_dict(five_element_pharmacophore())
    aromatic_2 = {}
    aromatic_2["name"] = "Aromatic"
    aromatic_2["hasvec"] = True
    aromatic_2["svector"] = {}
    aromatic_2["svector"]["x"] = 0.0
    aromatic_2["svector"]["y"] = 0.0
    aromatic_2["svector"]["z"] = 1.0
    aromatic_2["x"] = 0.0
    aromatic_2["y"] = 0.9999999999999999
    aromatic_2["z"] = 1.9999999999999998
    aromatic_2["radius"] = 0.9999999999999999
    aromatic_2["enabled"] = True
    aromatic_2["vector_on"] = 0
    aromatic_2["minsize"] = ""
    aromatic_2["maxsize"] = ""
    aromatic_2["selected"] = False

    donor = {}
    donor["name"] = "HydrogenDonor"
    donor["hasvec"] = True
    donor["svector"] = {}
    donor["svector"]["x"] = 0.0
    donor["svector"]["y"] = 0.7071067811865475
    donor["svector"]["z"] = 0.7071067811865475
    donor["x"] = 0.9999999999999999
    donor["y"] = 1.9999999999999998
    donor["z"] = 1.9999999999999998
    donor["radius"] = 0.9999999999999999
    donor["enabled"] = True
    donor["vector_on"] = 0
    donor["minsize"] = ""
    donor["maxsize"] = ""
    donor["selected"] = False

    hydrophobic = {}
    hydrophobic["name"] = "Hydrophobic"
    hydrophobic["hasvec"] = False
    hydrophobic["svector"] = {}
    hydrophobic["svector"]["x"] = 1
    hydrophobic["svector"]["y"] = 0
    hydrophobic["svector"]["z"] = 0
    hydrophobic["x"] = -0.9999999999999999
    hydrophobic["y"] = 1.9999999999999998
    hydrophobic["z"] = 1.9999999999999998
    hydrophobic["radius"] = 0.9999999999999999
    hydrophobic["enabled"] = True
    hydrophobic["vector_on"] = 0
    hydrophobic["minsize"] = ""
    hydrophobic["maxsize"] = ""
    hydrophobic["selected"] = False

    expected = {}
    expected["points"] = []
    expected["points"].append(acceptor)
    expected["points"].append(donor)
    expected["points"].append(hydrophobic)
    expected["points"].append(aromatic_1)
    expected["points"].append(aromatic_2)

    assert pharmer == expected
