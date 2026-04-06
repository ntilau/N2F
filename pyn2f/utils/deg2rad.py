"""Convert degrees to radians"""
import numpy as np


def deg2rad(deg):
    """
    Convert degrees to radians.
    
    Parameters
    ----------
    deg : float or array_like
        Angle in degrees
        
    Returns
    -------
    rad : float or ndarray
        Angle in radians
    """
    return (np.pi / 180) * deg
