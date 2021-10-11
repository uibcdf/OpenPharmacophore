import random
import string

def random_string(lenght):
    """ Generate a random string of lowercase letters and digits

        Parameters
        ----------
        lenght: int
            Lenght of the resulting string

        Returns
        -------
        str 
            The random string
    """
    if not isinstance(lenght, int):
        raise ValueError("Lenght must be an integer")
    characters = string.ascii_lowercase + string.digits
    result_str = "".join(random.choice(characters) for i in range(lenght))
    return result_str
