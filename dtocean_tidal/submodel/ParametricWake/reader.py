# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Thomas Roc, Mathew Topper
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
import platform

import numpy as np
import pandas as pd

#TR: integration of Sandia's model - beta 1.0
from read_db_mod import read_db

def read_database(data_path):
    
    """Converts the CDF database into many pandas table"""
    
    # Read and load CFD database into a dataframe at import level
    # Read Ct and TI info file
    try:
        
        # Allocate memory and build database
        reader = initiate_reader(data_path)
        database = create_database(reader, data_path)
        
    finally:
        
        # Deallocate memory
        delete_reader()
    
    return database
    
def create_database(reader, data_path):
    
    """Create pandas tables from database reader"""
    
    # Read and load CFD database into a dataframe at import level
    # Read Ct and TI info file
    cts = list(set([float(i) for i in open(os.path.join(data_path,
                                                        "Ct_set.txt"))]))
    cts.sort()
    ctsnames = ['ct'+str(i) for i in cts]

    tis = list(set([float(i) for i in open(os.path.join(data_path,
                                                        "TI_set.txt"))]))
    tis.sort()
    tisnames = ['ti'+str(i) for i in tis]
    
    #  Generate dataframe
    dfU = pd.DataFrame(index=ctsnames, columns=tisnames)
    dfV = pd.DataFrame(index=ctsnames, columns=tisnames)
    dfTKE = pd.DataFrame(index=ctsnames, columns=tisnames)
    
    for ii in tis:
        for cc in cts:
            reader.read_u_v_tke(cc, ii)
            dfU['ti'+str(ii)]['ct'+str(cc)] = np.copy(reader.u)
            dfV['ti'+str(ii)]['ct'+str(cc)] = np.copy(reader.v)
            dfTKE['ti'+str(ii)]['ct'+str(cc)] = np.copy(reader.tke)
            dfX = np.copy(reader.x)
            dfY = np.copy(reader.y)
    
    # Build the dictionaries        
    database = {'dfU'   : dfU,
                'dfV'   : dfV,
                'dfTKE' : dfTKE,
                'dfX'   : dfX,
                'dfY'   : dfY,
                'cts'   : np.asarray(cts),
                'tis'   : np.asarray(tis)}
    
    return database

def initiate_reader(data_path):

    """Bring up the reader module

    Returns:
      fortran object
    """

    # Set path to data folder
    # data_path = os.path.join(os.path.normpath(data_path), "")

    if not platform.system()=='Windows':
        read_db.datapath = "{:<512}".format(data_path+"/")
    else:
        read_db.datapath = "{:<512}".format(os.path.normpath(data_path) + "\\")

    # Allocate memory
    read_db.nx=800
    read_db.ny=320
    read_db.alloc_arrays()

    return read_db

def delete_reader():

    """Break down the reader module"""

    # Deallocate memory
    read_db.dealloc_arrays()

    return

