"""Plotting functions for near and far fields."""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ..utils import vector2matrix
from .plot_properties import get_figure_properties


def sf_plot_box_nf(matrix_size, nf_data, title='Box Near Field'):
    """
    Plots scalar near field on box surface.

    Parameters
    ----------
    matrix_size : tuple
        Shape of the field data matrix
    nf_data : ndarray
        Near field data to plot
    title : str, optional
        Plot title
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    X, Y = np.meshgrid(range(matrix_size[1]), range(matrix_size[0]))
    Z_mag = np.abs(nf_data.reshape(matrix_size))
    
    surf = ax.plot_surface(X, Y, Z_mag, cmap='viridis')
    ax.set_title(title)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('|NF|')
    fig.colorbar(surf)
    plt.show()


def sf_plot_sph_nf(matrix_size, nf_data, title='Spherical Near Field'):
    """
    Plots scalar near field on spherical surface.

    Parameters
    ----------
    matrix_size : tuple
        Shape of the field data matrix (theta, phi)
    nf_data : ndarray
        Near field data to plot
    title : str, optional
        Plot title
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    theta_idx, phi_idx = np.meshgrid(range(matrix_size[0]), range(matrix_size[1]))
    Z_mag = np.abs(nf_data.reshape(matrix_size))
    
    # Convert spherical to Cartesian for plotting
    R = Z_mag
    X = R * np.sin(theta_idx) * np.cos(phi_idx)
    Y = R * np.sin(theta_idx) * np.sin(phi_idx)
    Z = R * np.cos(theta_idx)
    
    surf = ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8)
    ax.set_title(title)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    fig.colorbar(surf)
    ax.set_box_aspect([1, 1, 1])
    plt.show()


def sf_plot_sph_nf_solid(matrix_size, nf_data, title='Spherical Near Field'):
    """
    Plots scalar near field on spherical surface with solid rendering.

    Parameters
    ----------
    matrix_size : tuple
        Shape of the field data matrix (theta, phi)
    nf_data : ndarray
        Near field data to plot
    title : str, optional
        Plot title
    """
    sf_plot_sph_nf(matrix_size, nf_data, title)


def sf_plot_ff_cut_planes(phi_cuts, theta, pattern, title='Far Field Cut Planes'):
    """
    Plots far field cut planes.

    Parameters
    ----------
    phi_cuts : ndarray
        Phi values for cut planes
    theta : ndarray
        Theta values [degrees]
    pattern : ndarray
        Field pattern magnitude for each (theta, phi_cut)
    title : str, optional
        Plot title
    """
    fig, axes = plt.subplots(1, len(phi_cuts), figsize=(5*len(phi_cuts), 4))
    
    if len(phi_cuts) == 1:
        axes = [axes]
    
    for idx, phi_val in enumerate(phi_cuts):
        ax = axes[idx]
        ax.plot(theta, 20 * np.log10(np.abs(pattern[:, idx]) + 1e-10))
        ax.set_xlabel('Theta [deg]')
        ax.set_ylabel('Magnitude [dB]')
        ax.set_title(f'φ = {phi_val}°')
        ax.grid(True)
    
    fig.suptitle(title)
    plt.tight_layout()
    plt.show()


def vf_plot_array_geom(array_pos):
    """
    Plots vector field array geometry.

    Parameters
    ----------
    array_pos : ndarray (3, N_a)
        Array element positions
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    ax.plot(array_pos[0, :], array_pos[1, :], array_pos[2, :], 
            '.k', markersize=20)
    ax.set_xlabel('x [λ]')
    ax.set_ylabel('y [λ]')
    ax.set_zlabel('z [λ]')
    ax.set_title('Vector Array Geometry')
    ax.set_box_aspect([1, 1, 1])
    plt.show()


def vf_plot_ff_3d(theta, phi, pattern, title='3D Far Field Pattern'):
    """
    Plots 3D far field pattern.

    Parameters
    ----------
    theta : ndarray
        Elevation angles [radians]
    phi : ndarray
        Azimuth angles [radians]
    pattern : ndarray (N_phi, N_theta)
        Field pattern magnitude
    title : str, optional
        Plot title
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    R = np.abs(pattern)
    X = R * np.sin(theta) * np.cos(phi)
    Y = R * np.sin(theta) * np.sin(phi)
    Z = R * np.cos(theta)
    
    surf = ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8)
    ax.set_title(title)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    fig.colorbar(surf)
    ax.set_box_aspect([1, 1, 1])
    plt.show()


def vf_plot_ff_cut_planes(phi_cuts, theta, et_pattern, ep_pattern, title='Vector Far Field Cut Planes'):
    """
    Plots vector far field cut planes.

    Parameters
    ----------
    phi_cuts : ndarray
        Phi values for cut planes
    theta : ndarray
        Theta values [degrees]
    et_pattern : ndarray
        Theta component field pattern
    ep_pattern : ndarray
        Phi component field pattern
    title : str, optional
        Plot title
    """
    fig, axes = plt.subplots(len(phi_cuts), 2, figsize=(10, 4*len(phi_cuts)))
    
    if len(phi_cuts) == 1:
        axes = axes.reshape(1, 2)
    
    for idx, phi_val in enumerate(phi_cuts):
        axes[idx, 0].plot(theta, 20 * np.log10(np.abs(et_pattern[:, idx]) + 1e-10))
        axes[idx, 0].set_ylabel('Magnitude [dB]')
        axes[idx, 0].set_title(f'Theta Pol., φ = {phi_val}°')
        axes[idx, 0].grid(True)
        
        axes[idx, 1].plot(theta, 20 * np.log10(np.abs(ep_pattern[:, idx]) + 1e-10))
        axes[idx, 1].set_ylabel('Magnitude [dB]')
        axes[idx, 1].set_title(f'Phi Pol., φ = {phi_val}°')
        axes[idx, 1].grid(True)
    
    axes[-1, 0].set_xlabel('Theta [deg]')
    axes[-1, 1].set_xlabel('Theta [deg]')
    fig.suptitle(title)
    plt.tight_layout()
    plt.show()


def vf_plot_ff_polar_cut_planes(phi_cuts, theta, pattern, title='Polar Far Field Pattern'):
    """
    Plots far field in polar coordinates.

    Parameters
    ----------
    phi_cuts : ndarray
        Phi values for cut planes
    theta : ndarray
        Theta values [radians]
    pattern : ndarray
        Field pattern magnitude
    title : str, optional
        Plot title
    """
    fig = plt.figure()
    
    for idx, phi_val in enumerate(phi_cuts):
        ax = fig.add_subplot(1, len(phi_cuts), idx+1, projection='polar')
        r = np.abs(pattern[:, idx])
        theta_deg = np.degrees(theta)
        ax.plot(theta_deg, 20 * np.log10(r + 1e-10))
        ax.set_title(f'φ = {phi_val}°')
        ax.grid(True)
    
    fig.suptitle(title)
    plt.tight_layout()
    plt.show()


def vf_plot_surface_power_density(matrix_size, power_density, title='Surface Power Density'):
    """
    Plots surface power density.

    Parameters
    ----------
    matrix_size : tuple
        Shape of power density matrix
    power_density : ndarray
        Power density on surface
    title : str, optional
        Plot title
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    X, Y = np.meshgrid(range(matrix_size[1]), range(matrix_size[0]))
    Z = power_density.reshape(matrix_size)
    
    surf = ax.plot_surface(X, Y, Z, cmap='hot')
    ax.set_title(title)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Power Density')
    fig.colorbar(surf)
    plt.show()


def plot_svd_error(errors, title='SVD Truncation Error'):
    """
    Plots SVD truncation error.

    Parameters
    ----------
    errors : ndarray
        Error values
    title : str, optional
        Plot title
    """
    fig, ax = plt.subplots()
    ax.semilogy(errors, 'o-')
    ax.set_xlabel('Coefficient Index')
    ax.set_ylabel('Error [V]')
    ax.set_title(title)
    ax.grid(True)
    plt.show()
