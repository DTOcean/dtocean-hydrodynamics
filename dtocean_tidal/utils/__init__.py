#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
from pkg_resources import get_distribution

#Local import
from distance_from_streamline import distance_from_streamline
from interpolation import interp_at_point

#Credentials
__version__ = get_distribution('dtocean-hydrodynamics').version
__authors__ = ['DTOcean Tidal Developers']
__licence__ = 'GNU Affero GPL v3.0'
__copyright__ = 'Copyright (c) 2014 DTOcean'

