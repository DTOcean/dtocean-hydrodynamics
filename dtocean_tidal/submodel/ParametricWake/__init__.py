#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
from pkg_resources import get_distribution

#Local import
from .wakeClass import WakeShape, Wake
from .reader import read_database

#Credentials
__version__ = get_distribution('dtocean-hydrodynamics').version
__authors__ = ['DTOcean Tidal Developers']
__licence__ = 'GNU Affero GPL v3.0'
__copyright__ = 'Copyright (c) 2014 DTOcean'

