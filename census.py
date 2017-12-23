import os
from collections import OrderedDict
import xml.etree.ElementTree as ElementTree
import json
import requests
import csv
import sqlite3

SCRIPT_DIR = os.path.dirname(__file__)
DATA_SUBDIR = os.path.join(SCRIPT_DIR, 'data', 'census')

with open(os.path.join(SCRIPT_DIR, 'census_fields.csv')) as f:
    CENSUS_FIELDS = [
        (fieldname, int(length_str), int(start_str) - 1)
        for fieldname, length_str, start_str in csv.reader(f)
    ]

def sf1_state_subdir(abbrev):
    return os.path.join(DATA_SUBDIR, 'sf1', abbrev)

def sf1_geo_state_db_filename(name, abbrev, fips_id):
    db_filename = os.path.join(sf1_state_subdir(abbrev), '{}geo2010_sf1.sqlite3'.format(abbrev))
    
    if not os.path.exists(db_filename):
        src_filename = sf1_geo_state_src_filename(name, abbrev, fips_id)
        assert os.path.exists(src_filename)
        
        with sqlite3.connect(db_filename) as db:
            db.execute('CREATE TABLE sf1 ({})'.format(
                ', '.join([field[0] for field in CENSUS_FIELDS])
            ))
            
            with open(src_filename) as f:
                for line in f:
                    vals = [
                        line[offset:offset + length]
                        for fieldname, length, offset in CENSUS_FIELDS
                    ]
                    db.execute('INSERT INTO sf1 VALUES ({})'.format(
                        ', '.join(['?'] * len(CENSUS_FIELDS))
                    ), vals)
    
    return db_filename

def sf1_geo_state_src_filename(name, abbrev, fips_id):
    src_filename = os.path.join(sf1_state_subdir(abbrev), '{}geo2010.sf1'.format(abbrev))
    assert os.path.exists(src_filename)
    
    # TODO: download file if not present

def vtd_state_shapefiles_filename(name, abbrev, fips_id):
    filename_base = 'tl_2012_{}_vtd10'.format(fips_id)
    dirname = os.path.join(DATA_SUBDIR, 'vtd', filename_base)
    filename = os.path.join(dirname, filename_base)
    assert os.path.isdir(dirname)
    assert os.path.exists('{}.shp'.format(filename))
    
    return filename
