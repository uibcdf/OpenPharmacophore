def load_smarts_fdef(file_name):
    """ Load custom feature definitions from a txt file. 

        Parameters
        ----------
        file_name : str
            Name of the file containing the smarts feature definitions
        
        Returns
        -------
        features : dict
            Dictionary which keys are SMARTS strings and values are feature names

        Notes
        -----
        The file must contain a SMARTS string follwed by the feature name. 
        
        Example:

        # Aromatic
        a1aaaaa1 Aromatic
        a1aaaa1 Aromatic

        Lines started with # are considered as comments
    """
    features = {} 
    # Load custom features file
    with open(file_name, "r") as file: 
        for line in file:
            if line[0] == "#" or line[0] == '\n':
                continue
            line = line.split(' ')
            feat_def = line[0]
            feat_name = line[1].rstrip()
            features[feat_def] = feat_name

    return features