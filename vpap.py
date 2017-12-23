import os

SCRIPT_DIR = os.path.dirname(__file__)
DATA_SUBDIR = os.path.join(SCRIPT_DIR, 'data', 'vpap')

def va_precinctmaps_filename(year):
    dirname = os.path.join(DATA_SUBDIR, 'va-precinct-maps-{}'.format(year), 'shp')
    assert os.path.isdir(dirname)
    filename = os.path.join(dirname, 'vaprecincts2016')
    assert os.path.exists('{}.shp'.format(filename))
    return filename

def va_manualprecinctlocations_filename(year):
    filename = os.path.join(DATA_SUBDIR, 'manual-precinct-locations-{}.sqlite'.format(year))
    assert os.path.exists(filename)
    return filename
