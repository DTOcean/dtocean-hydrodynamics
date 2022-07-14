# -*- coding: utf-8 -*-

#    Copyright (C) 2016-2022 Mathew Topper
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import os

import pytest

from polite.paths import Directory
from dtocean_hydro import start_logging
from dtocean_hydro.configure import get_install_paths


def test_start_logging():
    start_logging()


def test_get_install_paths_conda(mocker, tmpdir, install_lines):
    
    exe_path = tmpdir / "python.exe"
    ini_file = tmpdir / "etc" / "dtocean-data" / "install.ini"
    ini_file.write(install_lines, ensure=True)
    
    mocker.patch('polite.paths.sys.executable', new=str(exe_path))
    mocker.patch('polite.paths.system', new='win32')
    mocker.patch('dtocean_hydro.configure.SiteDataDirectory',
                 return_value=Directory(str(tmpdir)))
    
    paths = get_install_paths()
    
    assert set(['bin_path',
                'wec_share_path',
                'tidal_share_path']) == set(paths.keys())
    assert paths["bin_path"] == os.path.join("mock", "bin_mock")
    assert paths["wec_share_path"] == os.path.join("mock", "dtocean_wec_mock")
    assert paths["tidal_share_path"] == os.path.join("mock",
                                                     "dtocean_tidal_mock")


def test_get_install_paths_installer(mocker, tmpdir, install_lines):
    
    ini_file = tmpdir / "install.ini"
    ini_file.write(install_lines, ensure=True)
    
    mocker.patch('polite.paths.site_data_dir', return_value=str(tmpdir))
    mocker.patch('dtocean_hydro.configure.EtcDirectory',
                 return_value=Directory(str(tmpdir / "mock")))
    
    paths = get_install_paths()
    
    assert set(['bin_path',
                'wec_share_path',
                'tidal_share_path']) == set(paths.keys())
    assert paths["bin_path"] == os.path.join("mock", "bin_mock")
    assert paths["wec_share_path"] == os.path.join("mock", "dtocean_wec_mock")
    assert paths["tidal_share_path"] == os.path.join("mock",
                                                     "dtocean_tidal_mock")


def test_get_install_paths_missing(mocker, tmpdir):
    
    etc_path = tmpdir / "etc"
    site_path = tmpdir / "site"
    
    mocker.patch('dtocean_hydro.configure.EtcDirectory',
                 return_value=Directory(str(etc_path)))
    mocker.patch('dtocean_hydro.configure.SiteDataDirectory',
                 return_value=Directory(str(site_path)))
    
    with pytest.raises(RuntimeError) as excinfo:
        get_install_paths()
    
    assert str(etc_path) in str(excinfo)
    assert str(site_path) in str(excinfo)
