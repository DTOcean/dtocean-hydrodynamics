# -*- coding: utf-8 -*-


from __future__ import division

import os
import numpy as np
import pandas as pd
import argparse
import datetime as dt

from scipy.integrate import simps
from scipy.interpolate import interp1d
from shapely import wkb

def make_wave_statistics(wave_df,
                         period_bin_size=1.,
                         wave_height_bin_size=0.5,
                         direction_bin_size=45,
                         save_flag=False,
                         filepath='wave_stats'):
                             
    """
    Created on Sat Dec 20 09:45:30 2014
    Scatter diagram with directionality
    @author: fro
    
    IEC DTS 62600-101  IEC 2014
    pg - 33
    Te 1s
    Hm0 0.5 m
    dir 45 degrees

"""
                         
    dT = period_bin_size
    dH = wave_height_bin_size
    dD = direction_bin_size
    
    if (dT>1 or dH>0.5 or dD>45):
        
        errStr = ('ERROR: The given discretisation is larger then the'
                  'suggested by the IEC standard 101')
        raise ValueError(errStr)
        
    if 'Te' not in wave_df.columns:
        
        errStr = ('Enery period, Te, must be provided.')
        raise ValueError(errStr)
    
    # Calculate the number of bins in the direction space
    nDir = int(360. / dD)
    
    # Re-evaluate the number od direction in order to have a integer number of
    # bins
    dD = 360. / nDir
        
    wave_df = wave_df.dropna()
    
    direcs = (wave_df['Dir']) * np.pi / 180.
    Hm0_max = max(wave_df['Hm0'])
    Te_max = max(wave_df['Te'])
    nT = int(Te_max / dT)
    nH = int(Hm0_max / dH)
    
    # nT and nH are increased by two in order to have a full coverage of the
    # space.
    Tbin = np.array(range(nT + 2), dtype=float) * dT
    Hbin = np.array(range(nH + 2), dtype=float) * dH
    
    #discretisation in the angular dimension in degrees
    Dbin = list((np.array(range(nDir+1), dtype=float) / nDir * 360.
                                                            - (dD / 2.))
                                                              / 180. * np.pi)
    direcs[direcs>Dbin[-1]] -= 2 * np.pi
        
    thCut = pd.cut(direcs, Dbin, precision=6, include_lowest=True)
    
    wave_df['cuts'] = thCut[:] # thCut.labels
    
    dataGr = wave_df.groupby('cuts')
    
    SD = np.zeros((nH + 1, nT + 1, nDir))
    
    binns = [Hbin,Tbin]
    
    for ind, data2d in enumerate(dataGr):
        
        if data2d[1].empty:
            
            counts = 0
            
        else:

            D2d = np.array(zip(data2d[1]['Hm0'], data2d[1]['Te']))
                    
            (counts, xedges, yedges) = np.histogram2d(D2d[:,0],
                                                      D2d[:,1],
                                                      bins=binns)

        SD[:,:,ind] = counts
      
    SDn = SD / len(wave_df)
    
    Hm0_centered = yedges[:-1] + dH / 2.
    Te_centered = xedges[:-1] + dT / 2.
    Dir_centered = np.array(Dbin[0:-1], dtype=float) * 180. / np.pi + (dD / 2.)
    
    dictOut = {'Tp': Te_centered,
               'Hs': Hm0_centered,
               'B': Dir_centered,
               'p': SDn}
                    
    if save_flag:

        save_path = filepath + '_output.tsv'
        logmsg = 'Saving the statistical data to {}.'.format(save_path) 
        
        fid = open(save_path,'wb')
        fid.write('Meteocean Condition\n')
        fid.write(filepath + '\n')
        fid.write('\n')
        fid.write('Te range [s]\n')
        for te in Te_centered:
            fid.write('{}\t'.format(te))
        fid.write('\n')
        fid.write('Hm0 range [m]\n')
        for hm in Hm0_centered:
            fid.write('{}\t'.format(hm))
        fid.write('\n')
        fid.write('Dir range [degrees]\n')
        for di in Dir_centered:
            fid.write('{}\t'.format(di))
        fid.write('\n')
        fid.write('Porbability of occurence [h/h_year], Te-columns and '
                  'Hm0-rows\n')
        for di_in in range(nDir):
            fid.write('Angle\t{}\n'.format(Dir_centered[di_in]))
            for h_in in range(nH):
                for t_in in range(nT):
                    fid.write('{}\t'.format(SDn[h_in,t_in,di_in]))
                fid.write('\n')
  
        fid.close()
        
    return dictOut
    
def make_tide_statistics(dictinput):
    
    '''
     Function that selects the subset of current fields to be run, based on the
     statistical analysis described in Deliverable 2.4
     Written by JF Filipot and M Peray, October 1st, 2015
     contact: jean.francois.filipot@france-energies-marines.org
    
     function input:
     dictinput = { 'U':U, 'V':V, 'TI':TI, 'SSH':SSH, 't':t, 'xc':xc, 'yc':yc,
                   'x':x, 'y':y, 'ns':ns}
     U: east-west current fields (np.shape(U)=[nx,ny,nt])
     V: north-south current fields (np.shape(V)=[nx,ny,nt])
     x: east-west coordinate (len(x)=nx)
     y: north-south coordinate (len(y)=ny)
     t: time vector (len(t)=nt)
     xc: east-west coordinate of the location of interest (where the time
     series for the statistical analysis will be extracted, len(xc)=1).
     yc: north-east coordinate of the location of interest (where the time
     series for the statistical analysis will be extracted, len(xc)=1).
     ns number of scenarii to be played (len(ns)=1).
    
     function output:
    
     dictoutput = {'V': V, 'U':U, 'p':p, 'TI':TI, 'x':x, 'y':y, 'SSH':SSH}
     U: east-west current fields (np.shape(U)=[ny,nx,ns])
     V: north-south current fields (np.shape(V)=[ny,nx,ns])
     p: probability of occurence of each scenario (len(p)=ns)
     x: east-west coordinate (len(x)=nx)
     y: north-south coordinate (len(y)=ny)
     t: time vector (len(t)=ns)
    '''

    # extract data from dictionary
    uf=dictinput['U']
    vf=dictinput['V']
    TI=dictinput['TI']
    SSH=dictinput['SSH']
    t=dictinput['t']
    x=dictinput['x']
    y=dictinput['y']
    xc=dictinput['xc']
    yc=dictinput['yc']
    ns=dictinput['ns']
    
    # start by interpolation to get a time serie    
    u=uf[np.argmin(abs(x-xc)), np.argmin(abs(y-yc)), :]
    v=vf[np.argmin(abs(x-xc)), np.argmin(abs(y-yc)), :]

    # define increment and PDF axis
    du = 0.01
    dv = 0.01
    umax = np.max(u) + du
    umin = np.min(u) - du
    vmax = np.max(v) + dv
    vmin = np.min(v) - dv

    ltabu=np.ceil(np.array([umax-umin])/du) 
    ltabv=np.ceil(np.array([vmax-vmin])/dv)
    
    # then define each axis and center the bins

    tabu=np.linspace(umin, umax, ltabu)
    tabv=np.linspace(vmin, vmax, ltabv)
    sPDF = (ltabu-1, ltabv-1)
    PDF=np.zeros(sPDF)
    ub=np.zeros(ltabu-1)
    vb=np.zeros(ltabv-1)

    for iu in range(0,ltabu-1):
        for iv in range(0,ltabv-1):
            booluv=np.array([(u > tabu[iu]) &
                             (u <= tabu[iu+1]) & 
                             (v > tabv[iv]) & 
                             (v <= tabv[iv+1])])
            NPDF=booluv.astype(int)
            PDF[iu,iv]=np.sum((NPDF==1)*1)
            ub[iu]=(tabu[iu]+tabu[iu+1])/2.
            vb[iv]=(tabv[iv]+tabv[iv+1])/2.

    # PDF normalization

    PDF=PDF/len(u)/du/dv

    # Now find the direction with maximum variance

    um=np.zeros(ns)
    vm=np.zeros(ns)
    Pb=np.zeros(ns)

    if (np.max(ub) -np.min(ub)) > (np.max(vb) - np.min(vb)):

       tabc=np.linspace(np.min(ub), np.max(ub), ns)

    # Now get the bounds for the probability estimation (largest variance u)

       tabc1=np.zeros(ns)
       tabc2=np.zeros(ns)
       dc=tabc[1]-tabc[0]
       tabc1=tabc-dc/2
       tabc2=tabc+dc/2

       for ic in range(0,ns):
           imu=np.argmin(abs(tabc[ic]-tabu-1))
           # avoid division by 0
           if np.sum(PDF[imu,:]) == 0.:
               vm[ic] = 0.
           else:
               vm[ic]=np.sum(PDF[imu,:]*vb*dv)/np.sum(PDF[imu,:]*dv)
           um[ic]=tabu[imu]
           imu1=np.argmin(abs(tabc1[ic]-tabu))
           imu2=np.argmin(abs(tabc2[ic]-tabu))
           Pb[ic]=np.sum(
               np.sum(PDF[np.max((0,imu1)):np.min((ltabu,imu2)),:]*du)*dv)
   
    # Now get the bounds for the probability estimation (largest variance v)

    else:
        
       tabc=np.linspace(np.min(vb), np.max(vb), ns)

       tabc1=np.zeros(ns)
       tabc2=np.zeros(ns)
       dc=tabc[1]-tabc[0]
       tabc1=tabc-dc/2
       tabc2=tabc+dc/2

       for ic in range(0,ns):
           imv=np.argmin(abs(tabc[ic]-tabv-1))
           # avoid division by 0
           if np.sum(PDF[:,imv]) == 0.:
               um[ic] = 0.
           else:
               um[ic]=np.sum(PDF[:,imv]*ub*du)/np.sum(PDF[:,imv]*du)
           vm[ic]=tabv[imv]
           imv1=np.argmin(abs(tabc1[ic]-tabv))
           imv2=np.argmin(abs(tabc2[ic]-tabv))
           Pb[ic]=np.sum(
               np.sum(PDF[:,np.max((0,imv1)):np.min((ltabv,imv2))]*du)*dv)

    # Now get the time step with the closest characteristic
    # try to minimize the difference (um - u)+(vm-v)

    itmin=np.zeros(ns)
    for ic in range(0,ns):
        itmin[ic]=np.argmin(abs(um[ic]-u)+abs(vm[ic]-v))
 
    # convert itmin to integer
    
    ind=itmin.astype(int)   
    inds=ind.tolist()
    V=vf[:,:,inds]
    U=uf[:,:,inds]
    t=t[inds]
    p=Pb
    TI=TI[:,:,inds]
    SSH=SSH[:,:,inds]

    # remove 0 probability bins
    zero_pb = np.where(p == 0.)
    U = np.delete(U, zero_pb, 2)
    V = np.delete(V, zero_pb, 2)
    SSH = np.delete(SSH, zero_pb, 2)
    TI = np.delete(TI, zero_pb, 2)
    p = np.delete(p, zero_pb)
    t = np.delete(t, zero_pb)
    ns = ns - np.array(zero_pb).size

    # output
    dictoutput = {'U'    : U,
                  'V'    : V,
                  'SSH'  : SSH,
                  'TI'   : TI,
                  'x'    : x,
                  'y'    : y,
                  'p'    : p,
                  't'    : t,
                  'ns'   : ns
                  }
    
    return dictoutput
    
def make_power_histograms(device_power_pmfs,
                          bin_width=None):
    
    '''This function converts the hydrodynamics output into the array power
    output histogram required for Electrical analysis.
    
    '''
    
    if bin_width is None: bin_width = 0.1
        
    # Asses the maximum power in all the pmfs
    max_power = 0.
    
    for power_pmf in device_power_pmfs.itervalues():
        
        dev_max = power_pmf[:, 1].max()
        if dev_max > max_power: max_power = dev_max

    # Set the power bins to include the maximum power
    power_bins = np.arange(0, max_power + bin_width, bin_width)

    device_hists = {}
    
    for dev_id,  dev_power_pmf in device_power_pmfs.iteritems():
        
        hist, final_bins = np.histogram(dev_power_pmf[:,0], bins=power_bins)
        output_occurrence = sum_bins(hist, dev_power_pmf[:,1])
        bin_widths = [j-i for i, j in zip(power_bins[:-1], power_bins[1:])]
        device_hists[dev_id] = (output_occurrence / np.array(bin_widths),
                                final_bins)
        
    return device_hists
    
def add_Te(wave_df, gamma=3.3, drop_nan=True):
    
    given_cols = set(wave_df.columns)
    needed_cols = set(['Hm0', 'Tp'])
    
    if not needed_cols.issubset(given_cols):
        
        errStr = ("Columns names 'Hm0' and 'Tp' must be provided in given "
                  "dataframe.")
        raise ValueError(errStr)
    
    t_H = wave_df['Hm0']
    t_T = wave_df['Tp']
    
    t_Te = []

    for t_Hi, t_Ti in zip(t_H, t_T):
        
        t_S, t_w = make_JONSWAP(t_Hi, t_Ti, gamma)
        new_spectrum = make_spectra_analysis(t_S, t_w)        
        
        t_Te.append(new_spectrum['Tm_10'])
            
    new_df = wave_df.copy()
    new_df['Te'] = t_Te
    
    if drop_nan:
        
        new_df = new_df.dropna()
                
    return new_df
    
def make_JONSWAP(Hm0, Tp, gamma=3.3, w=-1, wc=-1):
    
    if wc<0:
        wc = 33. / Tp
    
    if w<0:
        w = np.linspace(0, wc, 257)
    
    g = 9.8063
    sa = 0.07
    sb = 0.09
    wp = 2. * np.pi / Tp
    
    # if sa and sb are different from 0.07 and 0.09 the scaling factor A
    # needs to be evaluated as:
    # >>A=(Hm0/g)**2/16/simps(S,w)
    A = 5.061 * Hm0**2 / Tp**4 * (1 - 0.287 * np.log(gamma))
    
    s = sb * np.ones(len(w))
    s[w<wp] = sa
       
    S = A * g**2 / w**5 * np.exp(-5. / 4. * (wp / w)**4) * \
                                gamma**np.exp(-0.5 * ((w / wp - 1.) / s)**2)
    S[0] = 0 
    
    return S, w
    
        
def make_spectra_analysis(S,w):

    g = 9.8063
    f = w / 2. / np.pi
    
    m0 = simps(S*(f)**0, w)
    m1 = simps(S*(f)**1, w)
    m2 = simps(S*(f)**2, w)
    m3 = simps(S*(f)**3, w)
    m4 = simps(S*(f)**4, w)
    
    m_1 = simps(S[f>0] / f[f>0], w[f>0])  # = m_1
    
    
    Hm0  = 4. * np.sqrt(m0) 
    Tm01 = m0 / m1
    Tm02 = np.sqrt(m0 / m2) 
    Tm24 = np.sqrt(m2 / m4)
    Tm_10 = m_1 / m0
    
    Tm12 = m1 / m2
    
    Tp   = 1. / f[S.argmax()]                             # peak period /length
    Ss   = 2. * np.pi * Hm0 / g / Tm02**2                 # Significant wave steepness
    Sp   = 2. * np.pi * Hm0 / g / Tp**2                   # Average wave steepness 
    Ka   = abs(simps(S * np.exp(1j * w * Tm02), w)) / m0  # groupiness factor
    
    # Quality control parameter 
    # critical value is approximately 0.02 for surface displacement records
    # If Rs>0.02 then there are something wrong with the lower frequency part 
    # of S.
    smooth = interp1d(f, S)    
    
    Rs = sum(smooth([0.0146 * 2. * np.pi,
                     0.0195 * 2. * np.pi,
                     0.0244 * 2. * np.pi])) / 3. / max(S)
                     
    #Second estimation of Tp    
    Tp2  = 2. * np.pi * simps(S**4, w) / simps(w * S**4, w)
    
    alpha1 = Tm24 / Tm02                  # m(3)/sqrt(m(1)*m(5))
    eps2   = np.sqrt(Tm01 / Tm12 - 1)          # sqrt(m(1)*m(3)/m(2)^2-1)
    eps4   = np.sqrt(1 - alpha1**2)          # sqrt(1-m(3)^2/m(1)/m(5))
    Qp     = 2. / m0**2 * simps(w * S**2, w)
    
    dic = {'m0':        m0,
           'm1':        m1,
           'm2':        m2,
           'm3':        m3,
           'm4':        m4,
           'm_1':       m_1,
           'Hm0':       Hm0,
           'Tm01':      Tm01,
           'Tm02':      Tm02,
           'Tm24':      Tm24,
           'Tm_10':     Tm_10,
           'Tm12':      Tm12,
           'Tp':        Tp,
           'Ss':        Ss,
           'Sp':        Sp,
           'Ka':        Ka,
           'Rs':        Rs,
           'Tp2':       Tp2,
           'alpha1':    alpha1,
           'eps2':      eps2,
           'eps4':      eps4,
           'Qp':        Qp}

    return dic
    
def add_Te_interface():
    
    '''Command line interface for add_Te.
    
    Example:
    
        To get help::
        
            $ add-Te -h
            
    '''
    
    epiStr = ('Francesco Ferri, Mathew Topper (c) 2015.')
              
    desStr = "Add Te time series to wave data containing Hm0 and Tp."

    parser = argparse.ArgumentParser(description=desStr,
                                     epilog=epiStr)
    
    parser.add_argument("path",
                        help=("path to file containing wave height and peak "
                              "period time series (excel or csv)"),
                        type=str)
                        
    parser.add_argument("-g", "--gamma",
                        help=("JONSWAP spectrum gamma parameter"),
                        type=float,
                        default=3.3)
                        
                                     
    args = parser.parse_args()
        
    file_path   = args.path
    gamma       = args.gamma
    
    # Build a data frame from the given file
    _, ext = os.path.splitext(file_path)
    
    if ext == ".csv":
        
        wave_df = pd.read_csv(file_path)
        
    elif ".xls" in ext:
        
        wave_df = pd.read_excel(file_path)
        
    else:
        
        errStr = "File must be either CSV or MS Excel format"
        raise ValueError(errStr)
        
    new_df = add_Te(wave_df, gamma)
    
    if ext == ".csv":
        
        new_df.to_csv(file_path, index=False)
        
    elif ".xls" in ext:
        
        new_df.to_excel(file_path, index=False)

    return

def sort_by_first(arrays):
    
    '''Sort the first array and apply same order to all arrays in list.
    
    args:
        arrays (list): List of numpy.NDArray objects
    
    returns:
        list of sorted arrays.
    
    '''

    order = np.argsort(arrays[0])
    result = [x[order] for x in arrays]

    return result

def sum_bins(bins, data):
    
    '''
    
    args:
        bins (array): number of items in each bin.
        data (array): data to be counted.
    
            
    returns:
        bin_sums
    
    '''

    start = 0
    bin_sums = []
    for items in bins:
        end = start + items
        bin_sums.append(data[start:end].sum())
        start = end
    
    return np.array(bin_sums)
    
def bathy_records_to_strata(bathy_records):
    
    """Convert the bathymetry layers table returned by the database into 
    Strata structure raw input"""
    
    bathy_table = pd.DataFrame.from_records(bathy_records, columns=[
                                                            "utm_point",
                                                            "depth",
                                                            "mannings_no",
                                                            "layer_order",
                                                            "sediment_type",
                                                            "initial_depth",
                                                            "total_depth",
                                                            "terminal_depth"])
    
    x = []
    y = []
    
    for point_hex in bathy_table["utm_point"]:
        point = wkb.loads(point_hex, hex=True)
        coords = list(point.coords)[0]
        x.append(coords[0])
        y.append(coords[1])
        
    bathy_table["x"] = x
    bathy_table["y"] = y
    
    layers = list(set(bathy_table["layer_order"]))
    layers.sort()
    
    depth_layers = []
    sediment_layers = []
    
    xi = list(set(x))
    yj = list(set(y))
    xi.sort()
    yj.sort()

    for z in layers:
        
        depths = []
        sediments = []
        
        for y in yj:
            
            d = []
            s = []
            
            for x in xi:
                
                point_df = bathy_table.loc[(bathy_table['x'] == x) &
                                           (bathy_table['y'] == y) &
                                           (bathy_table['layer_order'] == z)]
                
                if point_df.empty:
                    
                    d.append(None)
                    s.append(None)
                    
                else:
                    
                    d.append((point_df["depth"] -
                                        point_df["initial_depth"]).values[0])
                    s.append(point_df["sediment_type"].values[0])
                    
            depths.append(d)
            sediments.append(s)
            
        depth_layers.append(depths)
        sediment_layers.append(sediments)
        
    depth_array = np.swapaxes(np.array(depth_layers, dtype=float), 0, 2)
    sediment_array = np.swapaxes(np.array(sediment_layers), 0, 2)
    
    layer_names = ["layer {}".format(x) for x in layers]
    
    raw_strata = {"values": {"depth": depth_array,
                             "sediment": sediment_array},
                  "coords": [xi, yj, layer_names]}
    
    return raw_strata
    
def bathy_records_to_mannings(bathy_records):
    
    """Convert the bathymetry layers table returned by the database into 
    Strata structure raw input"""
    
    bathy_table = pd.DataFrame.from_records(bathy_records, columns=[
                                                            "utm_point",
                                                            "depth",
                                                            "mannings_no",
                                                            "layer_order",
                                                            "sediment_type",
                                                            "initial_depth",
                                                            "total_depth",
                                                            "terminal_depth"])
    
    x = []
    y = []
    
    for point_hex in bathy_table["utm_point"]:
        point = wkb.loads(point_hex, hex=True)
        coords = list(point.coords)[0]
        x.append(coords[0])
        y.append(coords[1])
        
    bathy_table["x"] = x
    bathy_table["y"] = y
    
    mannings_grid = []
    
    xi = list(set(x))
    yj = list(set(y))
    xi.sort()
    yj.sort()
        
    for y in yj:
        
        mannings_row = []
        
        for x in xi:
            
            point_df = bathy_table.loc[(bathy_table['x'] == x) &
                                       (bathy_table['y'] == y) &
                                       (bathy_table['layer_order'] == 1)]
            
            if point_df.empty:
                
                mannings_row.append(None)
                
            else:
                
                mannings_row.append(point_df["mannings_no"].values[0])
                
        mannings_grid.append(mannings_row)
        
    mannings_array = np.array(mannings_grid, dtype=float)
    
    mannings_raw = {"values": mannings_array,
                    "coords": [xi, yj]}
    
    return mannings_raw
    
def tidal_series_records_to_xset(tidal_records):
    
    """Convert the bathymetry layers table returned by the database into 
    Strata structure raw input"""
    
    tidal_table = pd.DataFrame.from_records(tidal_records, columns=[
                                                    'fk_farm_array',
                                                    'project_bathymetry_id',
                                                    'fk_farm_id',
                                                    'local_index',
                                                    'fk_site_id',
                                                    'utm_point',
                                                    'measure_date',
                                                    'measure_time',
                                                    'u',
                                                    'v',
                                                    'id',
                                                    'turbulence_intensity',
                                                    'ssh',
                                                    'fk_point_id'])
    
    tidal_table = point_to_xy(tidal_table)

    tidal_table["datetime"] = [dt.datetime.combine(date, time) for
                                date, time in zip(tidal_table["measure_date"],
                                                  tidal_table["measure_time"])]
                                              
    tidal_table = tidal_table.drop("measure_date", 1)
    tidal_table = tidal_table.drop("measure_time", 1)
    
    xi = list(set(tidal_table["x"]))
    yj = list(set(tidal_table["y"]))
    xi.sort()
    yj.sort()
    
    steps = list(set(tidal_table["datetime"]))
    steps.sort()

    u_steps = []
    v_steps = []
    ssh_steps = []
    ti_steps = []
    
    for t in steps:
        
        us = []
        vs = []
        sshs = []
        tis = []
        
        for y in yj:
            
            u = []
            v = []
            ssh = []
            ti = []
            
            for x in xi:
                
                point_df = tidal_table.loc[(tidal_table['x'] == x) &
                                           (tidal_table['y'] == y) &
                                           (tidal_table['datetime'] == t)]
                
                if point_df.empty:
                    u.append(None)
                    v.append(None)
                    ssh.append(None)
                    ti.append(None)
                else:
                    u.append(point_df["u"].values[0])
                    v.append(point_df["v"].values[0])
                    ssh.append(point_df["ssh"].values[0])
                    ti.append(point_df["turbulence_intensity"].values[0])
                    
            us.append(u)
            vs.append(v)
            sshs.append(ssh)
            tis.append(ti)
            
        u_steps.append(us)
        v_steps.append(vs)
        ssh_steps.append(sshs)
        ti_steps.append(tis)
            
    u_array = np.swapaxes(np.array(u_steps, dtype=float), 0, 2)
    v_array = np.swapaxes(np.array(v_steps, dtype=float), 0, 2)
    ssh_array = np.swapaxes(np.array(ssh_steps, dtype=float), 0, 2)
    ti_array = np.swapaxes(np.array(ti_steps, dtype=float), 0, 2)
        
    raw = {"values": {"U": u_array,
                      'V': v_array,
                      "SSH": ssh_array,
                      "TI": ti_array},
           "coords": [xi, yj, steps]}
    
    return raw

def point_to_xy(df, point_column="utm_point"):
    
    x = []
    y = []
    
    for point_hex in df[point_column]:
        point = wkb.loads(point_hex, hex=True)
        coords = list(point.coords)[0]
        x.append(coords[0])
        y.append(coords[1])
        
    df["x"] = x
    df["y"] = y
    
    return df

