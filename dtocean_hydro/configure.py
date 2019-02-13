# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Mathew Topper
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

# Standard Library
import logging

# Helpers for configuration files
from polite.paths import (ObjDirectory,
                          SiteDataDirectory,
                          UserDataDirectory,
                          DirectoryMap)
from polite.configuration import ReadINI

# Start logging
module_logger = logging.getLogger(__name__)


def get_install_paths():
    
    """Pick the necessary paths to configure the external files for the wave
    and tidal packages."""
    
    source_dir = ObjDirectory(__name__, "config")
    user_data = UserDataDirectory("dtocean_hydro", "DTOcean", "config")
    user_data_map = DirectoryMap(user_data, source_dir)
    
    install_src_name = "install.ini"
    
    # Check for bundled indicator file
    if source_dir.isfile(".bundled"):
        install_dst_name = "install_bundled.ini"
    else:
        install_dst_name = "install.ini"
    
    log_msg = ("Install configuration file name set to "
               "'{}'").format(install_dst_name)
    module_logger.debug(log_msg)
    
    user_data_map.safe_copy_file(install_src_name,
                                 "{}.txt".format(install_dst_name))
    user_ini_reader = ReadINI(user_data_map, install_dst_name)
    
    # Get the root path from the site data path.
    site_data = SiteDataDirectory("DTOcean Hydrodynamics", "DTOcean")
    site_ini_reader = ReadINI(site_data, install_dst_name)
    
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

