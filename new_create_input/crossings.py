"""
Computes the crossings and normals at crossings
between splines created from body points and grid lines
"""

import scipy.interpolate as interp
from input_utils import X_DIR, Y_DIR, Cross
import numpy as np


def _cross_to_namedtuple(entry):
    """ Function to improve readability """
    return Cross(entry[0], entry[1], entry[2], entry[3], entry[4])


def compute_crossings(list_of_bodies, xgrid, ygrid):
    """
    Computes all points where a spline obtained from the points in
    list_of_bodies crosses the grid formed by xgrid, ygrid.
    Returns a list of (x_crossing, y_crossing, dx_crossing, dy_crossing)
    containing all crossings for all bodies
    """
    crossings = []
    for body in list_of_bodies:
        crossings.extend(_compute_body_crossings(body, xgrid, ygrid))
    crossings = _remove_out_of_bounds(np.array(crossings), xgrid, ygrid)

    return [_cross_to_namedtuple(i) for i in crossings]


def _compute_body_crossings(body, xgrid, ygrid):
    """
    Computes all points where a spline obtained from the body points crosses
    the grid formed by xgrid,ygrid.
    Returns a list of (x_crossing, y_crossing, nx_crossing, ny_crossing)
    for one body
    """
    crossings = []

    spline, _ = interp.splprep(
        (body.xvals, body.yvals), s=0, per=body.is_closed
    )
    t_val, c_val, k_val = spline
    c_val_x, c_val_y = c_val
    for y in ygrid:
        shifted_spline = (t_val, (c_val_x, c_val_y-y), k_val)
        crossings.extend(
            _crossing_in_one_direction(spline, shifted_spline, X_DIR)
        )

    for x in xgrid:
        # Root finding only works in y direction
        # so we pass a reflected spline
        shifted_spline = (t_val, (c_val_y, c_val_x-x), k_val)
        crossings.extend(
            _crossing_in_one_direction(spline, shifted_spline, Y_DIR)
        )
    return crossings


def _crossing_in_one_direction(spline, shifted_spline, direction):
    """
    Computes where the shifted spline crosses zero
    and returns the corresponding position using the original spline
    """

    crossings = []
    roots = interp.sproot(shifted_spline)

    pos, der = [], []
    if(len(roots[1])):
        pos = np.array(interp.splev(roots[1], spline)).T
        der = np.array(interp.splev(roots[1], spline, der=1)).T
        der = [d/np.linalg.norm(d) for d in der]
    # normal = (dy, -dx)
    for p, d in zip(pos, der):
        crossings.append((p[0], p[1], d[1], -d[0], direction))
    return crossings


def _remove_out_of_bounds(crossings, xgrid, ygrid):
    xmin, xmax = np.min(xgrid), np.max(xgrid)
    ymin, ymax = np.min(ygrid), np.max(ygrid)
    dx = xgrid[1]-xgrid[0]
    dy = ygrid[1]-ygrid[0]
    tol = 1e-3
    crossings = list(filter(
        lambda c: (xmin-dx*tol < c[0] < xmax+dx*tol) and
        (ymin-dy*tol < c[1] < ymax+dy*tol), crossings)
    )
    return crossings


def _main():
    """ Basic usage example """
    import numpy as np
    import matplotlib.pyplot as plt
    import body_def

    t_vals = np.linspace(0, 2*np.pi, 200, endpoint=False)
    a_body = body_def.Body(0.5*np.cos(t_vals)+0.1,
                           0.5*np.sin(t_vals)+0.3, True)
    b_body = body_def.Body(0.5*np.cos(t_vals)+1.9,
                           0.5*np.sin(t_vals)+0.1, True)
    body_list = [a_body, b_body]

    xmin, xmax = (0, 2)
    ymin, ymax = (-1, 1)
    lines, cols = (50, 50)
    xgrid = np.linspace(xmin, xmax, cols)
    ygrid = np.linspace(ymin, ymax, lines)

    crossings = compute_crossings(body_list, xgrid, ygrid)
    p_crossings = np.array(crossings)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    X, Y = np.meshgrid(xgrid, ygrid)
    ax.scatter(X, Y, color='k', s=30, marker="+")
    for body in body_list:
        ax.plot(body.xvals, body.yvals)
    ax.scatter(p_crossings[:, 0], p_crossings[:, 1])
    ax.quiver(p_crossings[:, 0], p_crossings[:, 1], p_crossings[:, 2],
              p_crossings[:, 3], scale_units='xy', angles='xy', scale=10)
    ax.set_xlim(xmin-0.1, xmax+0.1)
    ax.set_ylim(ymin-0.1, ymax+0.1)
    ax.set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    _main()
