# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [3.2.0] - 2022-09-12

### Changed

-   Renamed dtocean_wec.submodule.utils.mesh_util module as 
    dtocean_wec.submodule.utils.mesh and refactored.

### Fixed

-   Fixed bugs that occur if opening wrong type of mesh file in dtocean-wec.

## [3.1.1] - 2022-07-15

### Fixed

-   Removed rogue print statement.

## [3.1.0] - 2022-07-14

### Changed

-   Changed the way that the non-code data is discovered. It now looks for
    the dtocean-data conda package or the installer retrieved from the
    https://github.com/DTOcean/dtocean-data repository.

### Removed

-   The paths to the non-code data can no longer be modified in the user
    configuration.

## [3.0.0] - 2021-10-12

### Added

-   Added loop to rotor interaction solver. This produces a more accurate
    result than the single iteration used previously.
-   Added check for additional grid orientations to find the maximum number of 
    devices that can be placed during an optimisation. One additional device is 
    added in each direction (in combination) and then assessed.

### Changed

-   Improved log messages
-   Now uses the dominant wake model for superposition of rotor wake
    velocity deficits. The same dominant wake is chosen for the turbulence
    kinetic energy factor calculation.
-   Replaced term "induction" to avoid confusion. Prefer "velocity" or
    "velocity coefficient".
-   Replaced Mannings number input with beta and alpha coefficients for use
    with Soulsby type velocity profile calculation.
-   Improved performance of rotor CFD data reading.
-   Renamed deg360_to_radpi function as bearing_to_radians.
-   Improved streamlines plot.

### Fixed

-   Fixed turbine angle of attack calculation.
-   Fixed minimum distance between turbines check.
-   Fixed rated power per rotor for MCT style device.
-   Fixed rotor streamline calculation travelling in wrong direction.
-   Fix bug in positioning of rotors for MCT style device.
-   Fixed expected device layout orientation in optimiser by making the 
    inter-column spacing refer to the x-direction and inter-row to the y. 
    Array rotation is then the angle of the oncoming flow minus 90 degrees.
-   Fixed optimiser start point estimator increasing scale when it should be 
    reducing it
-   Fix incorrect sorting of outputs for rotors by making names equal length
-   Fix bug which forced all rotors to share the same input data

## [2.0.0] - 2019-03-05

### Added

-   Added change log.
-   Added continuous integration configuration files.
-   Added output of calculated power matrix for the isolated device.
-   Instantaneous device power may no longer exceed the rated power.

### Changed

-   Changed wave modules mean power output per device to return actual power
    generated rather than the sum of all powers across the power matrix.
-   Optimised EnergyProduction function by reducing number of calls to
    block_diag. This has provided a one third reduction in run times.
-   Reduce memory consumption of wave calculations by using single precision 
    complex numpy arrays. This is necessary for solving large OEC farms with 8MB
    RAM.
-   A alternative configuration file is now used for the location of the
    hydrodynamic data files if the module is bundled into the installer or if
    installed from source.

### Removed

-   Removed pin of numpy at version 1.10.1 and updated bootstrap command in 
    setup.py, which now requires MinGW to be installed separately. Note that 
    PyQt and matplotlib must still be pinned due to incompatibility of later 
    versions with Python 2.
-   Removed installer code and data provided in the DTOcean data package.

### Fixed

-   Numerous PEP8 fixes.
-   Tidal module velocity profile generator switched to Manning's formulation as
    was incorrectly set to the Soulsby type.
-   Fixed bugs in the array layout optimiser code.
-   Fixed bug in calculation of device depths for wave module approximation
    test.
-   Fixed depreciation warning when sending arguments to setup.py test.
-   Refactored distance_from_streamline to improve readability and correct issue
    with streamlines travelling upstream rather than downstream.
-   Fixed issues with using non  -rectangular domain in the tidal module.
-   Fixed issues determining depth excluded zones with non  -rectangular
    domains.
-   NaNs are now set to zero in interp_at_point and edge cases are better
    handled.
-   Fixed confusing variable names for inputs of wave energy period and peak 
    period for wave energy calculations.
-   Fixed bug where the angle of attack for yawed tidal turbines was being
    incorrectly calculated.
-   Fixed tidal current streamline plotting (when debug flag is True).
-   Fixed bug in dtocean_waves' Directional class which reduced energy output
    from wave energy calculations.
-   Fixed bug in conversion of compass bearings to trig angles.
-   Fixed bug in array main direction setting in tidal energy calculations.


## [1.0.0] -   2017  -01  -05

### Added

-   Initial import of dtocean  -core from SETIS.
