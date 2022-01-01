def discretize(val, bins):
    """ Return the bin which the value belongs to.
    
        Parameters
        ----------
        val : float
            The vale to be binned
            
        bins : np.ndarray
            Array of bins. It has to be one dimensional and monotonic.
            
        Returns
        -------
        int
            The bin which the value belongs to.
    """
    
    for ii in range(bins.shape[0] - 1):
        if val == bins[ii]:
            return bins[ii]
        elif val > bins[ii] and val < bins[ii + 1]:
            if val > (bins[ii] + bins[ii + 1]) / 2:
                return bins[ii + 1]
            else:
                return bins[ii]
            
def double_bin_pharmacophore_graph(distance, bins, delta):
    """ Assign two bin values to the distance between pharmacophoric points.
    
        Parameters
        ----------
        distance : float
            The distance that will be binned.
            
        bins : np.ndarray
            Array of bins. It has to be one dimensional and monotonic.
         
        delta : float
            The tolerance from which a distance value is considered to belong to
            the lower and upper bin. It has to be a value between 0 and 0.5
            
        Returns
        -------
        2-tuple of int
            The two bins assigned to the distance.
    """
    
    for ii in range(bins.shape[0] - 1):
        if distance == bins[ii]:
            return (bins[ii], bins[ii])
        elif distance > bins[ii] and distance < bins[ii + 1]:
            if distance - bins[ii] > delta:
                return (bins[ii], bins[ii + 1])
            else:
                return (bins[ii], bins[ii])