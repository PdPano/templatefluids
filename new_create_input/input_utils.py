"""
Functions and definitions common to many parts of the implementation
"""

import collections
import numpy as np
Cross = collections.namedtuple('Cross', ["x", "y", "nx", "ny", "direction"])
ImagePoint = collections.namedtuple(
    'ImagePoint',
    ["image_x",
     "image_y",
     "x_border",
     "y_border",
     "nx_border",
     "ny_border"]
)
X_DIR = 0.0  # When discontinuity crosses y cte line
Y_DIR = 1.0  # When discontinuity crosses x cte line
INVALID = 100


def find_nearest_index(array, value):
    return (np.abs(array-value)).argmin()
