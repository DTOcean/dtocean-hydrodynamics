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

tidal_data = []

if not platform.system()=='Windows':
            
    tidal_data.append('submodel/ParametricWake/*.so')

else:

    tidal_data.append('submodel/ParametricWake/*.pyd')

wec_db_files = []    
for root, directories, filenames in os.walk("dtocean_wec\src"):
    for filename in filenames:
        file_path = os.path.join(root,filename)
        wec_db_files.append(file_path.replace("dtocean_wec\\", ""))
                                  
packageData={'dtocean_tidal': tidal_data,
             'dtocean_hydro': ['config/*.ini',
                               'config/*.yaml'],
             'dtocean_wec'  : wec_db_files}
                               
                               
class Bootstrap(Command):
    
    user_options = []

    def initialize_options(self):
        """Abstract method that is required to be overwritten"""

    def finalize_options(self):
        """Abstract method that is required to be overwritten"""

    def run(self):
    
        # Clean libraries
        clean = Cleaner("library files",
                        ['.pyd', '.so'])
        clean()
        
        # Setup paths
        start_dir = os.getcwd()
        temp_dir = tempfile.mkdtemp()
        mod_root = 'read_db_mod'
        dst_dir = os.path.join(local_dir,
                               "dtocean_tidal/submodel/ParametricWake")
        
        if not platform.system() == 'Windows':
            
            build_file_name = "{}.so".format(mod_root)
            src_file_path = os.path.join(local_dir, 'src', 'read_db.f90')
            # Build the file
            system_call = ("f2py -c {} -m {}").format(src_file_path,
                                                      mod_root)

        else:

            build_file_name = "{}.pyd".format(mod_root)
            src_file_path = os.path.join(local_dir, 'src', 'read_db.f90')
            # Build the file
            system_call = ("f2py -c {} -m {} "
                           "--compiler=mingw32").format(src_file_path,
                                                        mod_root)

        os.chdir(temp_dir)
        os.system(system_call)
        os.chdir(start_dir)

        # Move the file
        build_file_path = os.path.join(temp_dir, build_file_name)
        dst_file_path = os.path.join(dst_dir, build_file_name)
        shutil.move(build_file_path, dst_file_path)

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
      version='1.1.dev0',
      description='Hydrodynamics module for the DTOcean tools',
      author=('Francesco Ferri, '
              'Pau Mercade, '
              'Thomas Roc, '
              'Mathew Topper, '
              'Chris Chartrand'),
      author_email=('ff@civil.aau.dk, '
                    'pmr@civil.aau.dk, '
                    'thomas.roc@itpower.co.uk, '
                    'damm_horse@yahoo.co.uk, ' 
                    'ccchart@sandia.gov'),
      license = "GPLv3",
      setup_requires=['numpy', 'shapely'],
      packages=find_packages(),
      install_requires=['cma',
                        'descartes',
                        'h5py',
                        'matplotlib',
                        'numpy',
                        'pandas',
                        'polite>=0.9',
                        'pyopengl',
                        'scikit-learn',
                        'scipy',
                        'shapely',
#                        'PyQt4'
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
      
