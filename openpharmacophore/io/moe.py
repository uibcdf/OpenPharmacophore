from openpharmacophore import PharmacophoricPoint
from openpharmacophore.utils.bisection import insort_right
import pyunitwizard as puw
import datetime
import re
from typing import List

def from_moe(file_name: str) -> List[PharmacophoricPoint]:
    """ Loads a pharmacophore from a MOE ph4 file.

        Parameters
        ----------
        file_name : str
            Name of the file containing the pharmacophore.

        Returns
        -------
        points : list of openpharmacophore.PharmacophoricPoint
            A list of pharmacophoric points.
    """
    moe_to_oph = {
        "Cat": "positive charge",
        "Ani": "negative charge",
        "Don": "hb donor", 
        "Acc": "hb acceptor",
        "Hyd": "hydrophobicity",
        "Aro": "aromatic ring",
    }

    points = []

    with open(file_name, "r") as f:
        moe_ph4_string = f.read()
        
    # pieces = re.split("\s", moe_ph4_string)
    pieces = moe_ph4_string.split()
    pattern = r"Don|Acc|Aro|Cat|Ani|Hyd"
    # Index where the excluded volumes spheres coordinates, if any, start
    exclusion_inx = None
    for i, piece in enumerate(pieces):
        if re.match(pattern, piece):
            if len(piece) == 7: 
            # Double features have seven characters i.e Acc|Aro 
                feat_name_1 = piece[0: 3]
                feat_name_2 = piece[4:]
                # Coordinates are always two items after the feature name
                x = float(pieces[i + 2])
                y = float(pieces[i + 3])
                z = float(pieces[i + 4])
                radius = float(pieces[i + 5])
                point_1 = PharmacophoricPoint(
                    feat_type=moe_to_oph[feat_name_1],
                    center=puw.quantity([x, y, z], "angstroms"),
                    radius=puw.quantity(radius, "angstroms")
                )
                point_2 = PharmacophoricPoint(
                    feat_type=moe_to_oph[feat_name_2],
                    center=puw.quantity([x, y, z], "angstroms"),
                    radius=puw.quantity(radius, "angstroms")
                )
                insort_right(points, point_1, key=lambda p: p.short_name)
                insort_right(points, point_2, key=lambda p: p.short_name)
            else:
            # Regular features and
            # features with weird names i.e Cat$mDon, Acc2
                feat_name = piece[0:3]
                # Coordinates are always two items after the feature name
                x = float(pieces[i + 2])
                y = float(pieces[i + 3])
                z = float(pieces[i + 4])
                radius = float(pieces[i + 5])
                point = PharmacophoricPoint(
                    feat_type=moe_to_oph[feat_name],
                    center=puw.quantity([x, y, z], "angstroms"),
                    radius=puw.quantity(radius, "angstroms")
                ) 
                insort_right(points, point, key=lambda p: p.short_name)

        if piece == "#volumesphere":
            # Excluded volume spheres are 10 items after #volumesphere
            exclusion_inx = i + 10
            break
    
    if exclusion_inx:
        count = 1
        for p in pieces[exclusion_inx :]:
            if p == "#volume":
                break
            if count == 1:
                x = float(p)
                count += 1
            elif count == 2:
                y = float(p)
                count += 1
            elif count == 3:
                z = float(p)
                count += 1
            elif count == 4:
                radius = float(p)
                excluded_sphere = PharmacophoricPoint(
                    feat_type="excluded volume",
                    center=puw.quantity([x, y, z], "angstroms"),
                    radius=puw.quantity(radius, "angstroms")
                )
                insort_right(points, excluded_sphere, key=lambda p: p.short_name)
                count = 1
    
    return points
            
def _moe_ph4_string(pharmacophoric_points: List[PharmacophoricPoint]) -> str:
    """ Returns a string with the necessary info to create a MOE ph4 file.

        Parameters
        ----------
        pharmacophoric_points : List of PharmacophoricPoint
            List with the pharmacophoric points.

        file_name : str
            Name of the file containing the pharmacophore.

        Returns
        ------
        ph4_str : str
            The pharmacophore string.
    """
    oph_to_moe = {
        "aromatic ring": "Aro",
        "hydrophobicity": "Hyd",
        "hb acceptor": "Acc",
        "hb donor": "Don",
        "positive charge": "Cat",
        "negative charge": "Ani",
    }

    now = datetime.datetime.now()
    ph4_str = "#moe:ph4que" + " " + str(now.year) + "." + str(now.month) + "\n"
    ph4_str += "#pharmacophore 5 tag t value *\n"
    ph4_str += "scheme t Unified matchsize i 0 title t s $\n"
    ph4_str += f"#feature {len(pharmacophoric_points)} expr tt color ix x r y r z r r r ebits ix gbits ix\n"

    excluded_spheres = []
    for element in pharmacophoric_points:
        if element.feature_name == "excluded volume":
            excluded_spheres.append(element)
            continue
        feat_name = oph_to_moe[element.feature_name]
        center = puw.get_value(element.center, to_unit="angstroms")
        radius = str(puw.get_value(element.radius, to_unit="angstroms")) + " "
        x = str(center[0]) + " "
        y = str(center[1]) + " "
        z = str(center[2]) + " "
        ph4_str += feat_name + " df2f2 " + x + y + z + radius    
        ph4_str += "0 300 "
    
    if excluded_spheres:
        ph4_str += "\n#volumesphere 90 x r y r z r r r\n"
        for excluded in excluded_spheres:
            center = puw.get_value(excluded.center, to_unit="angstroms")
            radius = str(puw.get_value(excluded.radius, to_unit="angstroms")) + " "
            x = str(center[0]) + " "
            y = str(center[1]) + " "
            z = str(center[2]) + " "
            ph4_str += x + y + z + radius

    ph4_str += "\n#endpharmacophore"

    return ph4_str