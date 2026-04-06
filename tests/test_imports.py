"""
Test and verify Python module imports
"""

import sys
import os

# Add parent directory to path so n2f can be imported
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    # Test main module imports
    import pyn2f
    print("✓ Main pyn2f module imported successfully")
    
    # Test submodule imports
    from pyn2f import geometry, transforms, scalar, vector, utils, plotting
    print("✓ All submodules imported successfully")
    
    # Test key function imports from main module
    from pyn2f import (
        build_array, build_box, build_sphere,
        cartesian2spherical, spherical2cartesian,
        sf_nf_solver, vf_nf_solver,
        deg2rad, rad2deg,
        plot_sph_geom
    )
    print("✓ All key functions imported successfully from n2f")
    
    # Test numpy functionality
    import numpy as np
    arr = build_array(1, 3, 0.5, 3, 0.5)
    print(f"✓ build_array works: shape {arr.shape}")
    
    # Test deg2rad
    rad_val = deg2rad(90)
    expected_pi_2 = np.pi / 2
    assert np.abs(rad_val - expected_pi_2) < 1e-10
    print(f"✓ deg2rad(90) = {rad_val:.6f} rad (expected {expected_pi_2:.6f})")
    
    # Test rad2deg
    deg_val = rad2deg(np.pi / 2)
    assert np.abs(deg_val - 90) < 1e-10
    print(f"✓ rad2deg(π/2) = {deg_val:.6f}° (expected 90)")
    
    # Test vector2matrix
    test_vec = np.arange(12)
    mat = utils.vector2matrix((3, 4), test_vec)
    assert mat.shape == (3, 4)
    print(f"✓ vector2matrix works: (12,) -> {mat.shape}")
    
    print("\n✓✓✓ All tests passed! Library is working correctly. ✓✓✓")
    
except Exception as e:
    print(f"✗ Error during testing: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
