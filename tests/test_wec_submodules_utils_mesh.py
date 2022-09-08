
import os

import numpy as np
import pytest

from dtocean_wec.submodule.utils.mesh import (strip_comments,
                                              read_NEMOH,
                                              read_WAMIT)


@pytest.fixture
def cube_xsim():
    return 0


@pytest.fixture
def cube_vertices():
    return np.array([[ 1., -1.,  0.],
                     [ 1., -1., -2.],
                     [-1., -1., -2.],
                     [-1., -1.,  0.],
                     [-1., -1.,  0.],
                     [-1., -1., -2.],
                     [-1.,  1., -2.],
                     [-1.,  1.,  0.],
                     [-1.,  1.,  0.],
                     [-1.,  1., -2.],
                     [ 1.,  1., -2.],
                     [ 1.,  1.,  0.],
                     [ 1.,  1.,  0.],
                     [ 1.,  1., -2.],
                     [ 1., -1., -2.],
                     [ 1., -1.,  0.],
                     [ 1.,  1., -2.],
                     [-1.,  1., -2.],
                     [-1., -1., -2.],
                     [ 1., -1., -2.]])


def test_strip_comments():
    
    test = ("This is a multi-line string\n"
            "! With some whole line comments\n"
            "      ! And some blank space before a whole line comment\n"
            "And some comments that ! come after the token\n"
            "But there are not comments here.")
    expected = ("This is a multi-line string\n"
                "And some comments that\n"
                "But there are not comments here.")
    
    assert strip_comments(test) == expected


@pytest.mark.parametrize("file_name", ["cube.dat", "cube2.dat"])
def test_read_nemoh(test_data_folder,
                    file_name,
                    cube_xsim,
                    cube_vertices):
    
    cube_connectivity = np.array([[ 2,  3,  4,  1],
                                  [ 6,  7,  8,  5],
                                  [10, 11, 12,  9],
                                  [14, 15, 16, 13],
                                  [18, 19, 20, 17]], dtype=int) - 1
    
    file_path = os.path.join(test_data_folder, file_name)
    xsim, vertices, connectivity = read_NEMOH(file_path)
    
    assert xsim == cube_xsim
    assert (vertices == cube_vertices).all()
    assert (connectivity == cube_connectivity).all()


@pytest.mark.parametrize("file_name", ["cube.GDF", "cube2.GDF"])
def test_read_wamit(test_data_folder,
                    file_name,
                    cube_xsim,
                    cube_vertices):
    
    cube_connectivity = np.array([[ 0,  1,  2,  3],
                                  [ 4,  5,  6,  7],
                                  [ 8,  9, 10, 11],
                                  [12, 13, 14, 15],
                                  [16, 17, 18, 19]], dtype=int)
    
    file_path = os.path.join(test_data_folder, file_name)
    xsim, vertices, connectivity = read_WAMIT(file_path)
    
    assert xsim == cube_xsim
    assert (vertices == cube_vertices).all()
    assert (connectivity == cube_connectivity).all()
