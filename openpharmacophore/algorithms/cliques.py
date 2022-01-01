from openpharmacophore.pharmacophoric_point import distance_bewteen_pharmacophoric_points
from openpharmacophore.algorithms.discretize import discretize, double_bin_pharmacophore_graph
from openpharmacophore._private_tools.exceptions import OpenPharmacophoreValueError
import numpy as np
import pyunitwizard as puw

class PharmacophoreGraph():
    
    def __init__(self, pharmacophore, dmin=2.0, dmax=13.0, bin_size=1.0):
        
        self.adjacency_matrix = self.get_adjacency_matrix(pharmacophore, dmin, dmax, bin_size)
        
    @staticmethod
    def get_adjacency_matrix(pharmacophore, dmin=2.0, dmax=13.0, bin_size=1.0, delta=0.0):
        """ Compute the adjacency matrix of the pharmacophore.
        
            Parameters
            ----------
            pharmacophore : openpharmcophore.Pharmacophore
                A pharmacophore object.
            
            dmin : float
                The minimun distance in angstroms from which two pharmacophoric points are considered different.
            
            dmax : flaot
                The maximum distance in angstroms between pharmacohoric points.
                
            bin_size : float
                The size of the bins that will be used to bin the distances.
            
            delta : float
                The tolerance from which a distance value is considered to belong to
                the lower and upper bin. It has to be a value between 0 and 0.5
                
        """
        if dmin < 0 or dmax < 0:
            raise OpenPharmacophoreValueError("distance should be a positve number")
        if bin_size < 1:
            raise OpenPharmacophoreValueError("bin size must be a number greater or equal than 1")
        if delta < 0 or delta > 0.5:
            raise OpenPharmacophoreValueError("delta should be a number between 0 and 0.5")
        
        n_elements = pharmacophore.n_elements
        adj_matrix = np.empty((n_elements, n_elements), dtype=object)
        
        bins = np.arange(dmin, dmax, bin_size) 
       
        for ii in range(n_elements):
            for jj in range(ii, n_elements):
                if ii == jj:
                    adj_matrix[ii, jj] = pharmacophore.elements[jj].short_name
                else:
                    distance = distance_bewteen_pharmacophoric_points(
                            pharmacophore.elements[ii],
                            pharmacophore.elements[jj])
                    
                    if delta == 0.0:
                        binned_distance = discretize(distance, bins)
                        adj_matrix[ii, jj] = binned_distance
                        adj_matrix[jj, ii] = binned_distance
                    else:
                        bin_distance_1, bin_distance_2 = double_bin_pharmacophore_graph(distance, bins, delta)
                        adj_matrix[jj, ii] = bin_distance_1
                        adj_matrix[ii, jj] = bin_distance_2
                                            
        return adj_matrix
                
    def canonical_label(self):
        pass

def clique_detection_pharmacophore():
    pass