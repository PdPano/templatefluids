"""
Reads data from crossing to compute ghost point structures
"""

import numpy as np
import scipy.interpolate as interp
import scipy.spatial.distance as ssd
import scipy.optimize as optimize
from input_utils import X_DIR, Y_DIR, INVALID
import input_utils as iutils


def create_image_points(grid_mask, bodies, xgrid, ygrid):
    """
    Creates a dictionary of namedtuples containing:
    (gp_i, gp_j): (image_x, image_y, x_border, y_border, nx_border, ny_border)
    """
    image_point_pairs = {}
    splines = splines_from_bodies(bodies)
    if len(splines) > 0:
        i_inds, j_inds = np.where(grid_mask == -1)
        for i, j in zip(i_inds, j_inds):
            g_point = (xgrid[j], ygrid[i])
            spline, closest_u = _closest_spline(splines, g_point)
            image_point_pairs[(i, j)] = _create_image_point(
                g_point, spline, closest_u
            )
    return image_point_pairs


def _create_image_point(g_point, spline, closest_u):
    border_x, border_y = interp.splev(closest_u, spline)
    der_x, der_y = interp.splev(closest_u, spline, der=1)
    im_x, im_y = 2*border_x-g_point[0], 2*border_y-g_point[1]
    return iutils.ImagePoint(im_x[0], im_y[0],
                             border_x[0], border_y[0],
                             der_y[0], -der_x[0]
                             )


def _closest_spline(splines, g_point):
    """
    Find the spline closest to g_point and returns the spline and the spline
    parameter at the closest point
    """
    d_to_spline = np.inf
    closest_spline = None
    spline_u = -1
    for spline in splines:
        closest_u = optimize.fmin(
            _dist_to_p, x0=0.5, args=(g_point, spline), disp=0, xtol=1e-6
        )
        d_to_this_spline = _dist_to_p(closest_u, g_point, spline)
        if d_to_this_spline < d_to_spline:
            d_to_spline = d_to_this_spline
            closest_spline = spline
            spline_u = closest_u
    return closest_spline, spline_u


def create_grid_mask(crossings, xgrid, ygrid):
    """
    The resulting mask has:
      - -1 for ghost points
      - zero for fluid points
      - INVALID value for guard points inside the body (a big positive number)
    """
    grid_mask = np.zeros((len(ygrid), len(xgrid))).astype(int)
    for cross in crossings:
        if cross.direction == X_DIR:
            _safe_insert_x(grid_mask, cross, xgrid, ygrid)
        if cross.direction == Y_DIR:
            _safe_insert_y(grid_mask, cross, xgrid, ygrid)
    _fix_cross_der_points(grid_mask)
    _add_guards(grid_mask)
    grid_mask[grid_mask < 0] = -1
    grid_mask[(0 < grid_mask) & (grid_mask < INVALID)] = 0
    return grid_mask


def splines_from_bodies(bodies):
    splines = []
    for body in bodies:
        spline, _ = interp.splprep(
            (body.xvals, body.yvals), s=0, per=body.is_closed)
        splines.extend((spline,))
    return splines


def _dist_to_p(u, p, spline):
    s = interp.splev(u, spline)
    return ssd.euclidean(p, s)


def _safe_insert_x(grid_mask, cross, xgrid, ygrid):
    """Uses cross position to update grid_mask with bounds checking"""
    i, j = _get_index_x_crossing(cross, xgrid, ygrid)
    if j >= 0:
        grid_mask[i, j] += -1 if cross.nx >= 0 else 1
    if j+1 < len(xgrid):
        grid_mask[i, j+1] += 1 if cross.nx >= 0 else -1


def _safe_insert_y(grid_mask, cross, xgrid, ygrid):
    """Uses cross position to update grid_mask with bounds checking"""
    i, j = _get_index_y_crossing(cross, xgrid, ygrid)
    if i >= 0:
        grid_mask[i, j] += -1 if cross.ny >= 0 else 1
    if i+1 < len(ygrid):
        grid_mask[i+1, j] += 1 if cross.ny >= 0 else -1


def _get_index_x_crossing(cross, xgrid, ygrid):
    i = iutils.find_nearest_index(ygrid, cross.y)
    j = iutils.find_nearest_index(xgrid, cross.x)
    if cross.x-xgrid[j] < 0:
        j -= 1
    return i, j


def _get_index_y_crossing(cross, xgrid, ygrid):
    i = iutils.find_nearest_index(ygrid, cross.y)
    j = iutils.find_nearest_index(xgrid, cross.x)
    if cross.y-ygrid[i] < 0:
        i -= 1
    return i, j


def _fix_cross_der_points(grid_mask):
    imax, jmax = grid_mask.shape
    update_points = []
    diagonals = [(1, 1), (1, -1), (-1, 1), (-1, -1)]

    def check_bounds(i, j): return 0 <= i < imax and 0 <= j < jmax

    for i in range(imax):
        for j in range(jmax):
            for ishift, jshift in diagonals:
                if check_bounds(i+ishift, j+jshift):
                    if (grid_mask[i+ishift, j] < 0
                            and grid_mask[i, j+jshift] < 0
                            and grid_mask[i, j] == 0):
                        update_points.append((i, j))
    for p in update_points:
        grid_mask[p] = -1


def _add_guards(grid_mask):
    imax, jmax = grid_mask.shape
    update_points = []
    neighbors = [(1, 0), (-1, 0), (0, 1), (0, -1)]

    def check_bounds(i, j): return 0 <= i < imax and 0 <= j < jmax

    for i in range(imax):
        for j in range(jmax):
            if grid_mask[i, j] != 0:
                continue
            for ishift, jshift in neighbors:
                if check_bounds(i+ishift, j+jshift):
                    if grid_mask[i+ishift, j+jshift] < 0:
                        update_points.append((i, j))
    for p in update_points:
        grid_mask[p] = INVALID


def _create_plot_array_image_points(image_points, X, Y):
    x_plot = [None]
    y_plot = [None]
    for gpoint, image in image_points.items():
        x_plot.extend([X[gpoint], image.x_border, image.image_x, None])
        y_plot.extend([Y[gpoint], image.y_border, image.image_y, None])
    return x_plot, y_plot


def _main():
    """ Basic usage example """
    import crossings as cr
    import body_def
    import matplotlib.pyplot as plt
    # Creates two bodies
    t_vals = np.linspace(0, 2*np.pi, 200, endpoint=False)
    a_body = body_def.Body(0.5*np.cos(t_vals)-0.1,
                           0.5*np.sin(t_vals)+0.3, True)
    b_body = body_def.Body(0.5*np.cos(t_vals)+2.1,
                           0.5*np.sin(t_vals)+0.1, True)
    body_list = [a_body, b_body]

    # Grid definitions
    xmin, xmax = (0, 2)
    ymin, ymax = (-1, 1)
    lines, cols = (50, 50)
    xgrid = np.linspace(xmin, xmax, cols)
    ygrid = np.linspace(ymin, ymax, lines)
    X, Y = np.meshgrid(xgrid, ygrid)

    crossings = cr.compute_crossings(body_list, xgrid, ygrid)
    grid_mask = create_grid_mask(crossings, xgrid, ygrid)
    image_points = create_image_points(grid_mask, body_list, xgrid, ygrid)
    image_x_plot, image_y_plot = (
        _create_plot_array_image_points(image_points, X, Y)
    )

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(X, Y, color='k', s=30, marker="+")
    ax.scatter(X[grid_mask < 0], Y[grid_mask < 0], color='r')
    ax.scatter(X[grid_mask == INVALID], Y[grid_mask == INVALID], color='g')
    ax.scatter(np.array(crossings)[:, 0],
               np.array(crossings)[:, 1], marker="+")
    for body in body_list:
        ax.plot(body.xvals, body.yvals)
    ax.plot(image_x_plot, image_y_plot, marker='x')
    ax.set_xlim(xmin-0.1, xmax+0.1)
    ax.set_ylim(ymin-0.1, ymax+0.1)
    ax.set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    _main()
