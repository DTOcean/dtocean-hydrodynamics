# -*- coding: utf-8 -*-

#    Copyright (C) 2022 Mathew Topper
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

import os
import pytest

from polite.paths import Directory
from dtocean_wec.tab1 import ReadDb


@pytest.fixture
def readdb(mocker, qtbot, tmpdir, install_lines, main_window):
    
    exe_path = tmpdir / "python.exe"
    ini_file = tmpdir / "etc" / "dtocean-data" / "install.ini"
    ini_file.write(install_lines, ensure=True)
    
    mocker.patch('polite.paths.sys.executable', new=str(exe_path))
    mocker.patch('polite.paths.system', new='win32')
    mocker.patch('dtocean_hydro.configure.SiteDataDirectory',
                 return_value=Directory(str(tmpdir)))
    
    window = ReadDb(main_window)
    window.show()
    qtbot.addWidget(window)
    
    return window


def test_readdb_paths(qtbot, readdb):
    assert readdb._wec_db_folder == os.path.join("mock",
                                                 "dtocean_wec_mock",
                                                 "wec_db")
