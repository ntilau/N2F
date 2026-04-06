"""
Example: Near field to far field transformation from a bounding sphere
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '../../')

from pyn2f.geometry import build_array, get_sph_radius, build_sphere, geometry_utils
from pyn2f.utils import deg2rad
from pyn2f.scalar.sf_excitations import sf_excitations
from pyn2f.scalar import sf_solvers
from pyn2f.plotting import plot_sph_geom, sf_plot_ff_cut_planes

# Parameters for planar array of point sources on XY plane
lambda_ = 1
nbr_elems_x = 5  # number of point sources on x direction
wl_spacing_x = 0.5  # spacing between point sources in wavelengths (wl.) x dir
nbr_elems_y = 5  # analog on y
wl_spacing_y = 0.5  # analog on y

# Near field sphere parameters
flag = 0  # 0: dThetaNF-dPhiNF  1: sampling resolution 2: for nf plots
           # n.b. 2: induces errors in the n2f computation
sph_smpl_res = 0.05 * lambda_  # sampling resolution on the sphere
d_theta_nf = 5  # degrees
d_phi_nf = 5   # degrees
ext = 0.5  # minimum distance in wl. between the sphere and the array

# Building structure
print("Building array and sphere geometry...")
array_pos = build_array(lambda_, nbr_elems_x, wl_spacing_x,
                       nbr_elems_y, wl_spacing_y)

radius = get_sph_radius(lambda_, array_pos, ext)
sphere_pos, ds, theta_nf, phi_nf, matrix_size = build_sphere(
    radius, sph_smpl_res, d_theta_nf, d_phi_nf, flag)

r_mag, ndot_rv, n = geometry_utils.get_sph_vectors(array_pos, sphere_pos)

# Plot geometry
print("Plotting geometry...")
plot_sph_geom(matrix_size, array_pos, sphere_pos)

# Near field computation
print("Computing near field...")
steering_t = 45  # degrees
steering_p = 0   # degrees

excit_phasor = sf_excitations(lambda_, array_pos, steering_t, steering_p)
psi, del_psi = sf_solvers.sf_nf_solver(lambda_, excit_phasor, r_mag, ndot_rv)

# Near field to far field transformation
print("Computing far field...")
dtheta_ff = 0.5  # ff pattern resolution [degrees]
theta_ff = deg2rad(np.arange(-90, 91, dtheta_ff))
phi_ff = deg2rad([0, 90])

f_psi = sf_solvers.sf_nf2ff_solver(lambda_, theta_ff, phi_ff, sphere_pos, n, ds,
                                   psi, del_psi)
f_psi_ref = sf_solvers.sf_direct_ff_solver(lambda_, theta_ff, phi_ff,
                                          excit_phasor, array_pos)

# Compute gain
print("Computing gain...")
gain = sf_solvers.sf_compute_gain(f_psi)
ref_gain = sf_solvers.sf_compute_gain(f_psi_ref)

# Plot results
print("Plotting results...")
sf_plot_ff_cut_planes(theta_ff, gain, theta_ff, ref_gain, planes=[0, 1])
plt.show()

print("Done!")
