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
        elif bins[ii] < val < bins[ii + 1]:
            if val > (bins[ii] + bins[ii + 1]) / 2:
                return bins[ii + 1]
            else:
                return bins[ii]
