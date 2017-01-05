
from dtocean_hydro import start_logging
from dtocean_hydro.configure import get_install_paths

def test_get_install_paths():
    
    paths = get_install_paths()
    
    assert set(['bin',
                'wec_include',
                'tidal_include']) == set(paths.keys())
    assert isinstance(paths["bin"], str)
    
def test_start_logging():

    start_logging()
    
    assert True

