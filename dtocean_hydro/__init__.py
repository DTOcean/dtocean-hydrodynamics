#!/usr/bin/python2.7
# encoding: utf-8

import logging
from pkg_resources import get_distribution

from polite.paths import ObjDirectory, UserDataDirectory, DirectoryMap
from polite.configuration import Logger

# credentials
__version__ = get_distribution('dtocean-hydrodynamics').version
__authors__ = ['DTOcean Developers']

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
