# -*- coding: utf-8 -*-

#    Copyright (C) 2016-2018 Mathew Topper
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

from dtocean_hydro import start_logging
from dtocean_hydro.configure import get_install_paths

def test_get_install_paths():
    
    paths = get_install_paths()
    
    assert set(['bin',
                'wec_include',
                'tidal_include']) == set(paths.keys())
    assert isinstance(paths["bin"], str)
    
def test_start_logging():

    start_logging()
    
    assert True

