#!/usr/bin/python2.7
# encoding: utf-8
import os
import sys
import shutil
import tempfile
import platform

from setuptools import Command, find_packages, setup
from setuptools.command.test import test as TestCommand

local_dir = os.path.abspath(os.path.dirname(__file__))
mingw_dlls = ['libgcc_s_seh-1.dll',
              'libgfortran-4.dll',
              'libquadmath-0.dll',
              'libwinpthread-1.dll']
              
tidal_data = []

if not platform.system()=='Windows':
    tidal_data.append('submodel/ParametricWake/*.so')
else:
    tidal_data.append('submodel/ParametricWake/*.pyd')
    tidal_data.append('submodel/ParametricWake/*.dll')

packageData = {'dtocean_tidal': tidal_data,
               'dtocean_hydro': ['config/*.ini',
                                 'config/*.yaml',
                                 'config/.bundled']}


class Bootstrap(Command):
    
    user_options = []

    def initialize_options(self):
        """Abstract method that is required to be overwritten"""

    def finalize_options(self):
        """Abstract method that is required to be overwritten"""

    def run(self):
    
        # Clean libraries
        clean = Cleaner("library files",
                        ['.pyd', '.so', '.dll'])
        clean()
        
        # Setup paths
        start_dir = os.getcwd()
        temp_dir = tempfile.mkdtemp()
        mod_root = 'read_db_mod'
        
        src_file_path = os.path.join(local_dir, 'src', 'read_db.f90')
        dst_dir = os.path.join(local_dir,
                               "dtocean_tidal/submodel/ParametricWake")
        
        if not platform.system() == 'Windows':
            
            # Prepare build command
            build_file_name = "{}.so".format(mod_root)
            system_call = ("f2py -c {} -m {}").format(src_file_path,
                                                      mod_root)

        else:
        
            # Get mingw path from MINGW_BIN_PATH environment variable
            mingw_bin_path = os.getenv('MINGW_BIN_PATH')
            
            if mingw_bin_path is None:
                errStr = ("Environment variable MINGW_BIN_PATH must contain "
                          "path to folder containing mingw binaries")
                raise ValueError(errStr)

            # Add mingw to the path
            os.environ["PATH"] += os.pathsep + mingw_bin_path

            # Prepare build command
            build_file_name = "{}.pyd".format(mod_root)
            system_call = ("f2py -c {} -m {} "
                           "--compiler=mingw32").format(src_file_path,
                                                        mod_root)

        # Build the file
        os.chdir(temp_dir)
        os.system(system_call)
        os.chdir(start_dir)

        # Move the file
        build_file_path = os.path.join(temp_dir, build_file_name)
        dst_file_path = os.path.join(dst_dir, build_file_name)
        shutil.move(build_file_path, dst_file_path)
                
        # On Windows copy DLLs 
        if platform.system() == 'Windows': 
            for dll in mingw_dlls:
                dll_path = os.path.join(mingw_bin_path, dll)
                shutil.copy(dll_path, dst_dir)

        # Clean up the directory 
        os.removedirs(temp_dir)
        
        return
        
                     
class PyTest(TestCommand):

    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
    
        #import here, cause outside the eggs aren't loaded
        import pytest
        import shlex

        # Run the tests
        if self.pytest_args:
            opts = shlex.split(self.pytest_args)
        else:
            opts = []
        
        errno = pytest.main(opts)
        sys.exit(errno)


class CleanTest(Command):

    description = 'clean test files'
    user_options = []
     
    def initialize_options(self):
        pass
     
    def finalize_options(self):
        pass
     
    def run(self):
        clean = Cleaner("test files",
                        ['.pyc', '.pkl'])
        clean()
 
        
class Cleaner(object):
     
    def __init__(self, description='some files',
                       clean_list=None,
                       exclude_list=None
                       ):
        
        if clean_list is None: clean_list = []
        if exclude_list is None: exclude_list = ['.eggs',
                                                 '.git',
                                                 '.idea',
                                                 '.hg',
                                                 '__pycache__',
                                                 'test_data']
        
        self.description = description
        self.clean_list = clean_list
        self.exclude_list = exclude_list
     
    def is_exclude(self, path):
        for item in self.exclude_list:
            if path.find(item) != -1:
                return True
        return False
     
    def is_clean(self, path):
        return path.endswith(tuple(self.clean_list))
     
    def pickup_clean(self):
        for root, dirs, files in os.walk(os.getcwd()):
            if self.is_exclude(root):
                continue
            for fname in files:
                if not self.is_clean(fname):
                    continue
                yield os.path.join(root, fname)
                
    def __call__(self):
        print "start clean {}".format(self.description)
        for clean_path in self.pickup_clean():
            print "remove {}: {}".format(os.path.splitext(clean_path)[1],
                                         clean_path)
            os.remove(clean_path)
        print "end cleanup"


## Proceed to standard setup
setup(name='dtocean-hydrodynamics',
      version='2.0.0',
      description='Hydrodynamics module for the DTOcean tools',
      maintainer='Mathew Topper',
      maintainer_email='mathew.topper@dataonlygreater.com',
      license = "GPLv3",
      setup_requires=['numpy'],
      packages=find_packages(),
      install_requires=['cma',
                        'descartes',
                        'h5py',
                        'matplotlib<2',
                        'numpy',
                        'pandas',
                        'polite>=0.9',
                        'pyopengl',
#                        'PyQt4',
                        'scikit-learn',
                        'scipy',
                        'setuptools',
                        'shapely'
                        ],
      package_data=packageData,
      entry_points={
          'console_scripts':
              [
               'dtocean-wec = dtocean_wec:gui_interface',
               ]},
      zip_safe=False, # Important for reading config files
      tests_require=['pytest'],
      cmdclass={'bootstrap': Bootstrap,
                'test': PyTest,
                'cleantest': CleanTest}
      )
