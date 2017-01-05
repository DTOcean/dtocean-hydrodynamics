# -*- coding: utf-8 -*-
"""
"""

# Helpers for configuration files
from polite.paths import (ObjDirectory,
                          SiteDataDirectory,
                          UserDataDirectory,
                          DirectoryMap)
from polite.configuration import ReadINI

def get_install_paths(install_config_name="install.ini"):
    
    """Pick the necessary paths to configure the external files for the wave
    and tidal packages."""
    
    source_dir = ObjDirectory(__name__, "config")
    user_data = UserDataDirectory("dtocean_hydro", "DTOcean", "config")
    user_data_map = DirectoryMap(user_data, source_dir)
    user_data_map.safe_copy_file(install_config_name,
                                 "{}.txt".format(install_config_name))
    user_ini_reader = ReadINI(user_data_map, install_config_name)
    
    # Get the root path from the site data path.
    site_data = SiteDataDirectory("DTOcean Hydrodynamics", "DTOcean")
    site_ini_reader = ReadINI(site_data, install_config_name)
    
    if user_ini_reader.config_exists():
        config = user_ini_reader.get_config()
    elif site_ini_reader.config_exists():
        config = site_ini_reader.get_config()
    else:
        errStr = ("No suitable configuration file found at paths "
                  "{} or {}").format(site_ini_reader.get_config_path(),
                                     user_ini_reader.get_config_path())
        raise RuntimeError(errStr)

    path_dict = {"bin"              : config["dtocean_wec"]["bin_path"],
                 "wec_include"      : config["dtocean_wec"]["include_path"],
                 "tidal_include"    : config["dtocean_tidal"]["include_path"]
                 }

    return path_dict

