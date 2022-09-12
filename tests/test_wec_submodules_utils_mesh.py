
import os

import numpy as np
import pytest
import matplotlib.pyplot as plt

from dtocean_wec.submodule.utils.mesh import (strip_comments,
                                              read_NEMOH,
                                              read_WAMIT,
                                              MeshBem,
                                              _get_panel_norm,
                                              Panel)


@pytest.fixture
def cube_xsim():
    return 0


@pytest.fixture
def cube_vertices():
    return np.array([[ 1., -1.,  0.],
                     [-1., -1.,  0.],
                     [-1., -1., -2.],
                     [ 1., -1., -2.],
                     [-1., -1.,  0.],
                     [-1.,  1.,  0.],
                     [-1.,  1., -2.],
                     [-1., -1., -2.],
                     [-1.,  1.,  0.],
                     [ 1.,  1.,  0.],
                     [ 1.,  1., -2.],
                     [-1.,  1., -2.],
                     [ 1.,  1.,  0.],
                     [ 1., -1.,  0.],
                     [ 1., -1., -2.],
                     [ 1.,  1., -2.],
                     [ 1.,  1., -2.],
                     [ 1., -1., -2.],
                     [-1., -1., -2.],
                     [-1.,  1., -2.]])


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


@pytest.mark.parametrize("file_name", ["cube.GDF", "cube2.GDF"])
def test_mesh_bem_init_gdf(file_name,
                           test_data_folder,
                           cube_xsim,
                           cube_vertices):
    
    cube_connectivity = np.array([[ 0,  1,  2,  3],
                                  [ 4,  5,  6,  7],
                                  [ 8,  9, 10, 11],
                                  [12, 13, 14, 15],
                                  [16, 17, 18, 19]], dtype=int)
    mesh_bem = MeshBem(file_name, test_data_folder)
    
    assert mesh_bem.file_name == file_name
    assert mesh_bem.path == test_data_folder
    assert mesh_bem.xsim == cube_xsim
    assert (mesh_bem.vertices == cube_vertices).all()
    assert (mesh_bem.connectivity == cube_connectivity).all()
    assert mesh_bem.nP == len(cube_connectivity)
    assert mesh_bem.nV == len(cube_vertices)


@pytest.mark.parametrize("file_name", ["cube.dat", "cube2.dat"])
def test_mesh_bem_init_dat(file_name,
                           test_data_folder,
                           cube_xsim,
                           cube_vertices):
    
    cube_connectivity = np.array([[ 2,  3,  4,  1],
                                  [ 6,  7,  8,  5],
                                  [10, 11, 12,  9],
                                  [14, 15, 16, 13],
                                  [18, 19, 20, 17]], dtype=int) - 1
    mesh_bem = MeshBem(file_name, test_data_folder)
    
    assert mesh_bem.file_name == file_name
    assert mesh_bem.path == test_data_folder
    assert mesh_bem.xsim == cube_xsim
    assert (mesh_bem.vertices == cube_vertices).all()
    assert (mesh_bem.connectivity == cube_connectivity).all()
    assert mesh_bem.nP == len(cube_connectivity)
    assert mesh_bem.nV == len(cube_vertices)


def test_mesh_bem_bad_ext():
    
    with pytest.raises(IOError) as excinfo:
        MeshBem("bad.bad")
    
    assert "file type not supported" in str(excinfo.value)


@pytest.fixture
def mesh_bem(test_data_folder):
    return MeshBem("cube.GDF", test_data_folder)


@pytest.fixture
def mesh_bem_trans(mesh_bem):
    mesh_bem.translate(1, 1, -1)
    return mesh_bem


def test_mesh_bem_translate(mesh_bem_trans):
    assert (mesh_bem_trans.vertices == np.array([[ 2.,  0., -1.],
                                                 [ 0.,  0., -1.],
                                                 [ 0.,  0., -3.],
                                                 [ 2.,  0., -3.],
                                                 [ 0.,  0., -1.],
                                                 [ 0.,  2., -1.],
                                                 [ 0.,  2., -3.],
                                                 [ 0.,  0., -3.],
                                                 [ 0.,  2., -1.],
                                                 [ 2.,  2., -1.],
                                                 [ 2.,  2., -3.],
                                                 [ 0.,  2., -3.],
                                                 [ 2.,  2., -1.],
                                                 [ 2.,  0., -1.],
                                                 [ 2.,  0., -3.],
                                                 [ 2.,  2., -3.],
                                                 [ 2.,  2., -3.],
                                                 [ 2.,  0., -3.],
                                                 [ 0.,  0., -3.],
                                                 [ 0.,  2., -3.]])).all()


def test_mesh_bem_rotate(mesh_bem_trans):
    
    mesh_bem_trans.rotate(-np.pi / 2, (0, 0, 0))
    expected = np.array([[ 0., -2., -1.],
                         [ 0.,  0., -1.],
                         [ 0.,  0., -3.],
                         [ 0., -2., -3.],
                         [ 0.,  0., -1.],
                         [ 2.,  0., -1.],
                         [ 2.,  0., -3.],
                         [ 0.,  0., -3.],
                         [ 2.,  0., -1.],
                         [ 2., -2., -1.],
                         [ 2., -2., -3.],
                         [ 2.,  0., -3.],
                         [ 2., -2., -1.],
                         [ 0., -2., -1.],
                         [ 0., -2., -3.],
                         [ 2., -2., -3.],
                         [ 2., -2., -3.],
                         [ 0., -2., -3.],
                         [ 0.,  0., -3.],
                         [ 2.,  0., -3.]])
    
    assert np.isclose(mesh_bem_trans.vertices, expected).all()


def test_mesh_bem_invert_norm_all(mesh_bem):
    
    expected = np.array([[ 3,  2,  1,  0],
                         [ 7,  6,  5,  4],
                         [11, 10,  9,  8],
                         [15, 14, 13, 12],
                         [19, 18, 17, 16]])
    
    mesh_bem.invert_norm()
    assert (mesh_bem.connectivity == expected).all()


def test_mesh_bem_invert_norm_single(mesh_bem):
    
    expected = np.array([[ 3,  2,  1,  0],
                         [ 4,  5,  6,  7],
                         [ 8,  9, 10, 11],
                         [12, 13, 14, 15],
                         [16, 17, 18, 19]], dtype=int)
    
    mesh_bem.invert_norm(0)
    assert (mesh_bem.connectivity == expected).all()


def test_get_panel_norm_zpos():
    x = np.array([-1.,  1.,  1., -1.])
    y = np.array([-1., -1.,  1.,  1.])
    z = np.array([ 0.,  0.,  0.,  0.])
    norm = _get_panel_norm(x, y, z)
    assert np.isclose(norm, (0, 0, 1)).all()


def test_get_panel_norm_zneg():
    x = np.array([-1., -1.,  1., -1.])
    y = np.array([-1.,  1.,  1., -1.])
    z = np.array([ 0.,  0.,  0.,  0.])
    norm = _get_panel_norm(x, y, z)
    assert np.isclose(norm, (0, 0, -1)).all()


def test_get_panel_norm_yneg():
    x = np.array([ 1.,  1., -1., -1.])
    y = np.array([-1., -1., -1., -1.])
    z = np.array([-2.,  0.,  0., -2.])
    norm = _get_panel_norm(x, y, z)
    assert np.isclose(norm, (0, -1, 0)).all()


def test_panel():
    vertices = np.array([[-1., -1.,  0.],
                         [ 1., -1.,  0.],
                         [ 1.,  1.,  0.],
                         [-1.,  1.,  0.]])
    panel = Panel(vertices)
    assert np.isclose(panel.centroid, (0, 0, 0)).all()
    assert np.isclose(panel.n, (0, 0, 1)).all()


def test_panel_reverse():
    vertices = np.array([[-1., -1.,  0.],
                         [-1.,  1.,  0.],
                         [ 1.,  1.,  0.],
                         [ 1., -1.,  0.]])
    panel = Panel(vertices)
    assert np.isclose(panel.centroid, (0, 0, 0)).all()
    assert np.isclose(panel.n, (0, 0, -1)).all()


def test_panel_again():
    vertices = np.array([[ 1., -1., -2.],
                         [ 1., -1.,  0.],
                         [-1., -1.,  0.],
                         [-1., -1., -2.]])
    panel = Panel(vertices)
    assert np.isclose(panel.centroid, (0, -1, -1)).all()
    assert np.isclose(panel.n, (0, -1, 0)).all()


def test_mesh_bem_visualise_mesh(mocker, mesh_bem):
    mocker.patch("dtocean_wec.submodule.utils.mesh.plt.show")
    mesh_bem.visualise_mesh()
    plt.close("all")


def test_mesh_bem_visualise_norm(mocker, mesh_bem):
    mocker.patch("dtocean_wec.submodule.utils.mesh.plt.show")
    mesh_bem.visualise_norm()
    plt.close("all")


def test_mesh_bem_mesh_generation_gdf(tmpdir, test_data_folder, mesh_bem):
    
    mesh_bem.mesh_generation("gdf", str(tmpdir))
    files = tmpdir.listdir()
    
    assert len(files) == 1
    assert files[0].ext == ".gdf"
    
    expected = os.path.join(test_data_folder, "cube.GDF")
    
    with open(str(files[0])) as f1, open(expected) as f2:
        
        for _ in xrange(3):
            next(f1)
            next(f2)
        
        assert int(next(f1)) == int(next(f2))
        
        for line1, line2 in zip(f1, f2):
            vals1 = [float(x) for x in line1.split()]
            vals2 = [float(x) for x in line2.split()]
            assert np.isclose(vals1, vals2).all()


def test_mesh_bem_mesh_generation_nemoh(tmpdir, test_data_folder, mesh_bem):
    
    mesh_bem.mesh_generation("nemoh", str(tmpdir))
    files = tmpdir.listdir()
    
    assert len(files) == 1
    assert files[0].ext == ".dat"
    
    expected = os.path.join(test_data_folder, "cube.dat")
    
    with open(str(files[0])) as f1, open(expected) as f2:
        
        vals1 = [int(x) for x in next(f1).split()]
        vals2 = [int(x) for x in next(f2).split()]
        
        assert vals1 == vals2
        
        for line1, line2 in zip(f1, f2):
            
            vals1 = [float(x) for x in line1.split()]
            vals2 = [float(x) for x in line2.split()]
            
            if np.isclose(vals1, (0, 0, 0, 0)).all():
                break
            
            assert np.isclose(vals1, vals2).all()
        
        for line1, line2 in zip(f1, f2):
            vals1 = [int(x) for x in line1.split()]
            vals2 = sorted([int(x) for x in line2.split()])
            assert vals1 == vals2


def test_mesh_bem_mesh_generation_mesh(tmpdir, test_data_folder, mesh_bem):
    
    mesh_bem.mesh_generation("mesh", str(tmpdir))
    files = tmpdir.listdir()
    
    assert len(files) == 1
    assert files[0].ext == ".dat"
    
    expected = os.path.join(test_data_folder, "cube2.dat")
    
    with open(str(files[0])) as f1, open(expected) as f2:
        
        n_vertices = int(next(f1))
        
        assert n_vertices == int(next(f2))
        assert int(next(f1)) == int(next(f2))
        
        for _ in range(n_vertices):
            vals1 = [float(x) for x in next(f1).split()]
            vals2 = [float(x) for x in next(f2).split()]
            assert np.isclose(vals1, vals2).all()
        
        for line1, line2 in zip(f1, f2):
            vals1 = [int(x) for x in line1.split()]
            vals2 = sorted([int(x) for x in line2.split()])
            assert vals1 == vals2


def test_mesh_bem_mesh_generation_not_supported(mesh_bem):
    
    with pytest.raises(ValueError) as excinfo:
        mesh_bem.mesh_generation("bad")
    
    assert "Unsupported output format" in str(excinfo.value)
