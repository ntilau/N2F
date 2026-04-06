"""
Example: Near field to far field transformation from a bounding box
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '../../')

from pyn2f.geometry import build_array, get_box_dim, build_box, geometry_utils
from pyn2f.transforms import coordinate_transforms
from pyn2f.utils import deg2rad
from pyn2f.scalar.sf_excitations import sf_excitations
from pyn2f.scalar import sf_solvers
from pyn2f.plotting import sf_plot_ff_cut_planes

# Parameters
lambda_ = 1  # wavelength
nbr_elems_x = 3  # number of point sources on x direction
wl_spacing_x = 0.5  # spacing between point sources in wavelengths (x dir)
nbr_elems_y = 5  # number of point sources on y direction
wl_spacing_y = 0.5  # spacing between point sources in wavelengths (y dir)

# Near field box parameters
wl_ranging = 0.5
wl_spacing = 0.1
ext = 0.5

# Building structure
print("Building array and box geometry...")
array_pos = build_array(lambda_, nbr_elems_x, wl_spacing_x,
                       nbr_elems_y, wl_spacing_y)

(x_min, x_max, y_min, y_max, z_min, z_max, x_pts, y_pts, z_pts) = \
    get_box_dim(lambda_, array_pos, wl_ranging, wl_spacing, ext)

box_pos, box_n, ds, m_size = build_box(
    [1, 1, 1, 1, 1, 1], x_min, x_max, y_min, y_max, z_min, z_max,
    x_pts, y_pts, z_pts, 1, 0, 0)

r_mag, ndot_rv = geometry_utils.get_box_vectors(array_pos, box_pos, box_n)

# Steering parameters
steering_t = 0  # [degrees]
steering_p = 0  # [degrees]

print("Computing near field...")
excit_phasor = sf_excitations(lambda_, array_pos, steering_t, steering_p)
psi, del_psi = sf_solvers.sf_nf_solver(lambda_, excit_phasor, r_mag, ndot_rv)

# Near field to far field transformation
print("Computing far field...")
dtheta_ff = 1
theta_ff = deg2rad(np.arange(-90, 91, dtheta_ff))
phi_ff = deg2rad([0, 90])

f_psi = sf_solvers.sf_nf2ff_solver(lambda_, theta_ff, phi_ff, box_pos, box_n, ds,
                                   psi, del_psi)
f_psi_ref = sf_solvers.sf_direct_ff_solver(lambda_, theta_ff, phi_ff,
                                          excit_phasor, array_pos)

# Compute gain
gain = sf_solvers.sf_compute_gain(f_psi)
gain_ref = sf_solvers.sf_compute_gain(f_psi_ref)

# Plot results
print("Plotting results...")
sf_plot_ff_cut_planes(theta_ff, gain, theta_ff, gain_ref, planes=[0, 1])
plt.show()

print("Done!")
