# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [2.0.0] - 2019-03-05

### Added

- Added change log.
- Added continuous integration configuration files.
- Added output of calculated power matrix for the isolated device.
- Instantaneous device power may no longer exceed the rated power.

### Changed

- Changed wave modules mean power output per device to return actual power
  generated rather than the sum of all powers across the power matrix.
- Optimised EnergyProduction function by reducing number of calls to
  block_diag. This has provided a one third reduction in run times.
- Reduce memory consumption of wave calculations by using single precision 
  complex numpy arrays. This is necessary for solving large OEC farms with 8MB
  RAM.
- A alternative configuration file is now used for the location of the
  hydrodynamic data files if the module is bundled into the installer or if
  installed from source.

### Removed

- Removed pin of numpy at version 1.10.1 and updated bootstrap command in 
  setup.py, which now requires MinGW to be installed separately. Note that PyQt 
  and matplotlib must still be pinned due to incompatibility of later versions 
  with Python 2.
- Removed installer code and data provided in the DTOcean data package.

### Fixed

- Numerous PEP8 fixes.
- Tidal module velocity profile generator switched to Manning's formulation as
  was incorrectly set to the Soulsby type.
- Fixed bugs in the array layout optimiser code.
- Fixed bug in calculation of device depths for wave module approximation test.
- Fixed depreciation warning when sending arguments to setup.py test.
- Refactored distance_from_streamline to improve readability and correct issue
  with streamlines travelling upstream rather than downstream.
- Fixed issues with using non-rectangular domain in the tidal module.
- Fixed issues determining depth excluded zones with non-rectangular domains.
- NaNs are now set to zero in interp_at_point and edge cases are better
  handled.
- Fixed confusing variable names for inputs of wave energy period and peak 
  period for wave energy calculations.
- Fixed bug where the angle of attack for yawed tidal turbines was being
  incorrectly calculated.
- Fixed tidal current streamline plotting (when debug flag is True).
- Fixed bug in dtocean_waves' Directional class which reduced energy output from
  wave energy calculations.
- Fixed bug in conversion of compass bearings to trig angles.
- Fixed bug in array main direction setting in tidal energy calculations.


## [1.0.0] - 2017-01-05

### Added

- Initial import of dtocean-core from SETIS.
