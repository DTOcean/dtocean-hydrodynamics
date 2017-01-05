# DTOcean Hydrodynamics Module

## Installation

### Install Data Package

The hydrodynamic data package can be downloaded from
[Dropbox](https://www.dropbox.com/s/7i1sgp5j2y5vhog/dtocean-hydrodynamics-data-0.3.exe?dl=0).

Once downloaded execute the file to install. Remember to uninstall any old
versions first using the uninstaller in the DTOcean Hydrodynamics start menu
program folder.

### Set up an Anaconda environment

Using a windows command prompt enter the following commands:

```
conda create -n dtocean python pip pytest ipython-notebook
```

then

```
activate dtocean
```

or

```
C:\Anaconda\Scripts\activate.bat dtocean
```

### Add public Anaconda Cloud channel

```
conda config --add channels https://conda.anaconda.org/topper
```

### Install the dependencies

```
conda install configobj cma descartes h5py libpython=1.0 matplotlib numpy=1.10.1 pandas pyopengl pyqt=4.11.4 pyyaml scikit-learn scipy shapely-win-py27
```

The package "polite" must currently be installed separately from the "dev"
brach of the hg repository https://bitbucket.org/team_dtocean/polite using the
commands

```
cd path\to\polite
winmake.bat install
```

### Install the hydrodynamics package

```
cd path\to\dtocean-hydrodynamics
winmake.bat bootstrap
```

### Testing

```
cd ..\examples
python Ex_tidal.py
```

## DTOcean Tidal Array Submodule

### Introduction

The main aim of this tool is to take the hydrodynamic scenarios forward
to the succeeding WPs through resource assessment and definition of the
array layouts and the parameters in play for the array optimization. It
focuses on the hydrodynamic interaction between the individual devices
in the array, how this affects the resource, power performance, cost
uncertainties and environmental impact.

### Contacts

- Thomas Roc: thomas.roc@itpower.co.uk
- Vincent Neary: vsneary@sandia.gov
- Chris Chartrand: ccchart@sandia.gov
- Budi Gunawan: bgunawa@sandia.gov
- Jesse Roberts: jdrober@sandia.gov
- Jean-Francois Filipot: jean.francois.filipot@france-energies-
  marines.org

## DTOcean Project

DTOcean - "Optimal Design Tools for Ocean Energy Arrays" is funded by the 
European Commission’s 7th Framework Programme. Grant agreement number: 608597
