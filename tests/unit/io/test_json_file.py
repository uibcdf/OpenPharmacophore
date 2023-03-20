from openpharmacophore.io.json_file import json_pharmacophoric_elements, load_json_pharmacophore
import numpy as np
import pyunitwizard as puw
import math


def test_load_json_pharmacophore(json_pharmacophore_path):
    points, molecular_system, ligand = load_json_pharmacophore(
        json_pharmacophore_path, load_mol_sys=False)

    assert len(points) == 5
    assert molecular_system is None
    assert ligand is None

    assert points[0].feature_name == "hb acceptor"
    assert points[0].has_direction
    assert puw.get_value(points[0].radius, "angstroms") == 1.0
    assert np.allclose(puw.get_value(points[0].center, "angstroms"),
                       np.array([21.352, -14.531, 19.625]))
    assert np.allclose(points[0].direction,
                       np.array([-0.6405836470264256, 0.7029084735090229, -0.3091476492414897]))

    assert points[1].feature_name == "hb acceptor"
    assert points[1].has_direction
    assert puw.get_value(points[1].radius, "angstroms") == 1.0
    assert np.allclose(puw.get_value(points[1].center, "angstroms"),
                       np.array([19.355, -18.32, 23.987]))
    assert np.allclose(points[1].direction,
                       np.array([0.6859059711903811, 0.09092493673854565, 0.721987295292979]))

    assert points[2].feature_name == "hb donor"
    assert points[2].has_direction
    assert puw.get_value(points[2].radius, "angstroms") == 1.0
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
    assert puw.get_value(points[4].radius, "angstroms") == 2.0
    assert np.allclose(puw.get_value(points[4].center, "angstroms"),
                       np.array([19.985, -19.404402, 22.8422]))


def create_json_pharmacophoric_element(name, coords, radius, direction=None):

    has_direction = direction is not None
    if has_direction:
        direction = {
            "x": direction[0],
            "y": direction[1],
            "z": direction[2],
        }
    else:
        direction = {
            "x": 1.,
            "y": 0.,
            "z": 0.,
        }

    return {
        "name": name,
        "hasvec": has_direction,
        "svector": direction,
        "x": coords[0],
        "y": coords[1],
        "z": coords[2],
        "radius": radius,
        "enabled": True,
        "vector_on": 0,
        "minsize": "",
        "maxsize": "",
        "selected": False
    }


def test_json_pharmacophoric_elements(three_element_pharmacophore):

    # Expected output from to_json
    expected = {"points": []}

    aromatic = create_json_pharmacophoric_element(
        name="Aromatic",
        coords=(1.0, 0.0, 0.0),
        radius=1.0,
        direction=(0.0, 0.0, 1.0)
    )
    acceptor = create_json_pharmacophoric_element(
        name="HydrogenAcceptor",
        coords=(1.0, 2.0, 2.0),
        radius=1.0,
        direction=(0.0, 1 / math.sqrt(2.0), 1 / math.sqrt(2.0))
    )
    excluded = create_json_pharmacophoric_element(
        name="ExclusionSphere",
        coords=(2.0, 1.0, 2.0),
        radius=1.0
    )

    expected["points"].append(acceptor)
    expected["points"].append(excluded)
    expected["points"].append(aromatic)

    assert json_pharmacophoric_elements(three_element_pharmacophore) == expected
