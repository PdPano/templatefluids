"""
Defines basic body container
"""

import copy


class Body:
    def __init__(self, xvals, yvals, is_closed):
        self.xvals = copy.copy(xvals)
        self.yvals = copy.copy(yvals)
        self.is_closed = copy.copy(is_closed)
