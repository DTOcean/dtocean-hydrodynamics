# -*- coding: utf-8 -*-

#    Copyright (C) 2016-2022 Mathew Topper
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

import os
import logging

from polite.paths import (EtcDirectory,
                          SiteDataDirectory)
from polite.configuration import ReadINI

module_logger = logging.getLogger(__name__)


def get_install_paths():
    """Pick the necessary paths to configure the external files for the wave
    and tidal packages."""
    
    # Look in the etc directory
    etc_data = EtcDirectory("dtocean-data")
    etc_ini_reader = ReadINI(etc_data, "install.ini")
    
    # Get the root path from the possible site data paths
    site_data = SiteDataDirectory("DTOcean Data", "DTOcean")
    site_ini_reader = ReadINI(site_data, "install.ini")
    
    if etc_ini_reader.config_exists():
        config = etc_ini_reader.get_config()
    elif site_ini_reader.config_exists():
        config = site_ini_reader.get_config()
    else:
        errStr = ("No suitable configuration file found at paths "
                  "{} or {}").format(etc_ini_reader.get_config_path(),
                                     site_ini_reader.get_config_path())
        raise RuntimeError(errStr)
    
    prefix = config["global"]["prefix"]
    bin_path = os.path.join(prefix, config["global"]["bin_path"])
    wec_share_path = os.path.join(prefix, config["dtocean_wec"]["share_path"])
    tidal_share_path = os.path.join(prefix,
                                    config["dtocean_tidal"]["share_path"])
    
    return {"bin_path" : bin_path,
            "wec_share_path" : wec_share_path,
            "tidal_share_path" : tidal_share_path}
