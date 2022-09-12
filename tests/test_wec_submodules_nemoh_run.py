
import numpy as np
import pytest
import logging

from dtocean_wec.submodule.nemoh_run import _get_cylinder_radius
from dtocean_wec.submodule.utils.mesh import MeshBem


def test_get_cylinder_radius_no_meshs():
    
    with pytest.raises(RuntimeError) as excinfo:
        _get_cylinder_radius([])
    
    assert "Can not determine mesh extents" in str(excinfo)


def test_get_cylinder_radius_cube(test_data_folder):
    
    mesh = MeshBem("cube.GDF", test_data_folder)
    test = _get_cylinder_radius([mesh])
    expected = 2 + np.sqrt(2)
    
    assert np.isclose(test, expected)


def test_get_cylinder_radius_cube_offset(caplog, test_data_folder):
    
    mesh = MeshBem("cube.GDF", test_data_folder)
    mesh.translate(1, 0, 0)
    
    with caplog.at_level(logging.WARNING):
        test = _get_cylinder_radius([mesh])
    
    expected_d = 2 + np.sqrt(5)
    expected_log = "centroid of the mesh file is not centered"
    
    assert np.isclose(test, expected_d)
    assert expected_log in caplog.record_tuples[0][2]


def test_get_cylinder_radius_two_cubes(caplog, test_data_folder):
    
    mesh1 = MeshBem("cube.GDF", test_data_folder)
    mesh1.translate(1, 0, 0)
    
    mesh2 = MeshBem("cube.GDF", test_data_folder)
    mesh2.translate(-1, 0, 0)
    
    with caplog.at_level(logging.WARNING):
        test = _get_cylinder_radius([mesh1, mesh2])
    
    expected_d = 2 + np.sqrt(5)
    
    assert np.isclose(test, expected_d)
    assert not caplog.record_tuples
