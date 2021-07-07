# -*- coding: utf-8 -*-

#    Copyright (C) 2016-2021 Mathew Topper
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

import logging
from pkg_resources import get_distribution

from polite.paths import ObjDirectory, UserDataDirectory, DirectoryMap
from polite.configuration import Logger

# credentials
__authors__ = ['DTOcean Developers']
__version__ = get_distribution('dtocean-hydrodynamics').version

# Set default logging handler to avoid "No handler found" warnings.
try:  # Python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass

logging.getLogger(__name__).addHandler(NullHandler())

def start_logging(level=None):

    """Start python logger"""

    objdir = ObjDirectory(__name__, "config")
    datadir = UserDataDirectory("dtocean_hydro", "DTOcean", "config")
    dirmap = DirectoryMap(datadir, objdir)

    log = Logger(dirmap)
    log("dtocean_hydro",
        level=level,
        info_message="Begin logging for dtocean_hydro.")
