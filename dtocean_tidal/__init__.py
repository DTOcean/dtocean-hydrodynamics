#!/usr/bin/python2.7
# encoding: utf-8

import logging
from pkg_resources import get_distribution

# Local import
from dtocean_tidal.submodel.ParametricWake import read_database

# credentials
__version__ = get_distribution('dtocean-hydrodynamics').version
__authors__ = ['DTOcean Tidal Developers']
__licence__ = 'GNU Affero GPL v3.0'
__copyright__ = 'Copyright (c) 2014 DTOcean'

# Set default logging handler to avoid "No handler found" warnings.
try:  # Python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass

logging.getLogger(__name__).addHandler(NullHandler())

