"""
Handles the bodies and boundaries
"""
import numpy as np
from body_def import Body
from boundary_def import Boundary
import ast


class Geometry:
    def __init__(self):
        self.boundaries = []
        self.bodies = []
        self.lines = None
        self.columns = None
        self.xlims = None
        self.ylims = None

    def add_geometry(self, key, value):
        if key == "GEO_LINES":
            self.lines = ast.literal_eval(value)
            if not isinstance(self.lines, int):
                raise TypeError("Lines must be integer")
        if key == "GEO_COLUMNS":
            self.columns = ast.literal_eval(value)
            if not isinstance(self.lines, int):
                raise TypeError("Columns must be integer")
        if key == "GEO_XLIMS":
            xlims = ast.literal_eval(value)
        if key == "GEO_YLIMS":
            self.ylims = ast.literal_eval(value)

    def add_boundary(self, new_boundary):
        in_boundary = self._read_boundary(new_boundary)
        self.boundaries.append(in_boundary)

    def add_body_from_file(self, new_body):
        new_body_file, is_closed = new_body.split(',')
        new_body_file, is_closed = new_body_file.strip(), is_closed.strip()
        is_closed = {'TRUE': True, 'FALSE': False}[is_closed]
        body_vals = np.loadtxt(new_body_file, delimiter=',', skiprows=1)
        self.bodies.append(Body(body_vals[:, 0], body_vals[:, 1], is_closed))

    def _read_boundary(self, new_boundary):
        return Boundary(ast.literal_eval(new_boundary))
