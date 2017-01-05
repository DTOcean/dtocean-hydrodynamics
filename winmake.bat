@ECHO OFF
SET curdir=%cd%

IF "%1" == "bootstrap" CALL :Bootstrap
IF "%1" == "install" CALL :Install
IF "%1" == "test" CALL :Test

:Return
CD %curdir%
EXIT /B %ERRORLEVEL%

:: Function to bootstrap the source code
:Bootstrap
CD %~dp0 && (
    CALL python setup.py bootstrap
    CALL :Install
) || ( GOTO Return )
EXIT /B 0

:: Function to install the program
:Install
CD %~dp0 && (
    CALL python setup.py clean --all
    CALL python setup.py install
    CALL :Test
) || ( GOTO Return )
EXIT /B 0

:: Function to test the program
:Test
CD %~dp0 && (
    CALL python setup.py cleantest
    CALL python setup.py test -a "-rsx tests"
) || ( GOTO Return )
EXIT /B 0
