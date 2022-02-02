from openpharmacophore import PharmacophoricPoint
from openpharmacophore._private_tools.exceptions import InvalidFeatureError
import bisect

class PharmacophoricPointContainer():
    """ Class to store pharmacopohric points and keep them sorted by
        feature type. 
    """

    def __init__(self) -> None:
        self._points = []
        self._n_points = 0

    def add_point(self, pharmacophoric_point : PharmacophoricPoint) -> None:
        """ Add a pharmacophoric point to the container, keeping them sorted by feature type"""
        bisect.insort(self._points, pharmacophoric_point, key=lambda p : p.short_name)
    
    def remove_point(self, index : int) -> None:
        self._points.pop(index)
        self.n_pharmacophoric_points -=1

    def remove_points(self, indices : list) -> None:
        new_pharmacophoric_points = [element for i, element in enumerate(self._points) if i not in indices]
        self._points = new_pharmacophoric_points
        self.n_pharmacophoric_points = len(self._points)

    def remove_feature(self, feat_type : str) -> None:

        feats = PharmacophoricPoint.get_valid_features()
        if feat_type not in feats:
            raise InvalidFeatureError(f"Cannot remove feature. \"{feat_type}\" is not a valid feature type")

        tmp_points = [element for element in self._points if element.feature_name != feat_type]
        if len(tmp_points) == self.n_pharmacophoric_points: # No element was removed
            raise InvalidFeatureError(f"Cannot remove feature. The pharmacophore does not contain any {feat_type}")
        self._points = tmp_points
        self.n_pharmacophoric_points = len(self._points)

    def get_number_points(self) -> int:
        return self._n_points

    def tolist(self) -> list:
        return self._points

    def count_features(self) -> dict:
        counter = {
            "aromatic ring":   0,
            "hydrophobicity":  0,
            "hb acceptor":     0,
            "hb donor":        0,
            "positive charge": 0,
            "negative charge": 0,
        }

        for point in self._points:
            counter[point.feature_name] += 1

    def clear(self) -> None:
        self._points.clear()
        self._n_points = 0

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(n_points : {self._n_points})"