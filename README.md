[![appveyor](https://ci.appveyor.com/api/projects/status/github/DTOcean/dtocean-hydrodynamics?branch=master&svg=true)](https://ci.appveyor.com/project/DTOcean/dtocean-hydrodynamics)
[![codecov](https://codecov.io/gh/DTOcean/dtocean-hydrodynamics/branch/master/graph/badge.svg)](https://codecov.io/gh/DTOcean/dtocean-hydrodynamics)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/bb34506cc82f4df883178a6e64619eaf)](https://www.codacy.com/project/H0R5E/dtocean-hydrodynamics/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=DTOcean/dtocean-hydrodynamics&amp;utm_campaign=Badge_Grade_Dashboard&amp;branchId=8410911)
[![release](https://img.shields.io/github/release/DTOcean/dtocean-hydrodynamics.svg)](https://github.com/DTOcean/dtocean-hydrodynamics/releases/latest)

# DTOcean Hydrodynamics Module

This package provides the Hydrodynamics design module for the DTOcean tools.
It can calculate the energy output of arrays of fixed or floating wave or tidal
ocean energy converters, including the effect of interactions. It can optimise
the position of the devices for maximum energy yield, constrained by the given 
environment.

See [dtocean-app](https://github.com/DTOcean/dtocean-app) or [dtocean-core](
https://github.com/DTOcean/dtocean-app) to use this package within the DTOcean
ecosystem.

* For python 2.7 only.

## Installation

Installation and development of dtocean-hydrodynamics uses the [Anaconda 
Distribution](https://www.anaconda.com/distribution/) (Python 2.7)

These installation instructions are for WINDOWS ONLY.

### Data Package (Required for ALL installation Methods)

The hydrodynamic data package must be installed and can be downloaded from
[Github](https://setis.ec.europa.eu/dt-ocean/).

Once downloaded execute the file to install. If upgrading from version 1,
uninstall the old version first from the Windows start menu program folder,
or using the control panel. For version 2 and beyond, the uninstaller will
automatically remove the older version.

### Conda Package

To install:

```
$ conda install -c dataonlygreater dtocean-hydrodynamics
```

### Source Code

Conda can be used to install dependencies into a dedicated environment from
the source code root directory:

```
conda create -n _dtocean_hydro python=2.7 pip
```

Activate the environment, then copy the `.condrc` file to store installation  
channels:

```
$ conda activate _dtocean_hydro
$ copy .condarc %CONDA_PREFIX%
```

Install [polite](https://github.com/DTOcean/polite) into the environment. For 
example, if installing it from source:

```
$ cd \\path\\to\\polite
$ conda install --file requirements-conda-dev.txt
$ pip install -e .
```

Some modules in dtocean-hydrodynamics must be compiled before installation.
These instructions use [mingw-w64](https://mingw-w64.org), which should be
installed first. When ready, set the MINGW_BIN_PATH environment variable to
the "bin" folder of the mingw-w64 installation. For example:

```
$ SET MINGW_BIN_PATH=C:\mingw\mingw64\bin
```

Check that the variable is set correctly:

```
$ echo %MINGW_BIN_PATH%
C:\mingw\mingw64\bin
```

Now compile the modules:

```
$ cd \\path\\to\\dtocean-hydrodynamics
$ python setup.py bootstrap
```

Finally, install dtocean-hydrodynamics and its dependencies using conda and pip:

```
$ conda install --file requirements-conda-dev.txt
$ pip install -e .
```

To deactivate the conda environment:

```
$ conda deactivate
```

### Tests

A test suite is provided with the source code that uses [pytest](
https://docs.pytest.org).

If not already active, activate the conda environment set up in the [Source 
Code](#source-code) section:

```
$ conda activate _dtocean_hydro
```

Install pytest to the environment (one time only):

```
$ conda install -y pytest
```

Run the tests:

``` 
$ py.test tests
```

### Uninstall

To uninstall the conda package:

```
$ conda remove dtocean-hydrodynamics
```

To uninstall the source code and its conda environment:

```
$ conda remove --name _dtocean-hydro --all
```

To uninstall the data package use the link in the Windows start menu program
folder, or use the control panel.

## Usage

Example scripts are available in the "examples" folder of the source code.

For tidal energy converters:

```
cd examples
python Ex_tidal_v2.py
```

For wave energy converters:

```
cd examples
python Ex_wave_v2.py
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first to
discuss what you would like to change.

See [this blog post](
https://www.dataonlygreater.com/latest/professional/2017/03/09/dtocean-development-change-management/)
for information regarding development of the DTOcean ecosystem.

Please make sure to update tests as appropriate.

## Credits

This package was initially created as part of the [EU DTOcean project](
https://www.dtoceanplus.eu/About-DTOceanPlus/History) by:

 * Mathew Topper at [TECNALIA](https://www.tecnalia.com)
 * Vincenzo Nava at [TECNALIA](https://www.tecnalia.com)
 * Adam Colin at [the University of Edinburgh](https://www.ed.ac.uk/)
 * David Bould at [the University of Edinburgh](https://www.ed.ac.uk/)
 * Rui Duarte at [France Energies Marines](https://www.france-energies-marines.org/)
 * Francesco Ferri at [Aalborg University](https://www.en.aau.dk/)

It is now maintained by Mathew Topper at [Data Only Greater](
https://www.dataonlygreater.com/).

## License

[GPL-3.0](https://choosealicense.com/licenses/gpl-3.0/)
