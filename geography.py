import os
import sys
import shapefile
import sqlite3
import census
import geometry
import vpap
from districting import *

SCRIPT_DIR = os.path.dirname(__file__)
sys.path.append(SCRIPT_DIR)

def iter_precinct_shaperecords():
    sf = shapefile.Reader(census.vtd_state_shapefiles_filename('Virginia', 'va', '51'))
    for sr in sf.shapeRecords():
        yield sr.shape, {field[0] : value for field, value in zip(sf.fields[1:], sr.record)}
    sf.close()

def iter_precincts_2010():
    vtd_pop_dict = load_vtd_populations()
    
    for i, (shape, record) in enumerate(iter_precinct_shaperecords()):
        polygon = list(shape.points)
        if polygon[0] != polygon[-1]:
            polygon.append(polygon[0])
        
        yield Region(
            i, polygon,
            population = vtd_pop_dict[(record['STATEFP10'], record['COUNTYFP10'], int(record['VTDST10']))],
            record = record
        )

def load_precinct_map(precincts2010):
    # Because as of now it's annoying to get accurate maps and populations for 2017
    # precincts, we map 2017 precincts onto 2010 precincts and form districts using those.
    index = RegionIndex(precincts2010)
    precinct_map = {}
    #centroids = {}
    
    # Map 2016 precinct centroids to 2010 precincts
    sf = shapefile.Reader(vpap.va_precinctmaps_filename(2016))
    for sr in sf.shapeRecords():
        record = {field[0] : value for field, value in zip(sf.fields[1:], sr.record)}
        id = record['ID']
        polygon = [espg3857_to_wgs84(*point) for point in sr.shape.points]
        if polygon[0] != polygon[-1]:
            polygon.append(polygon[0])
        centroid = polygon_centroid(polygon)
        #centroids[id] = centroid
        precinct = index.region_for_point(centroid)
        precinct_map[id] = precinct
    sf.close()
    
    # Map manually added 2017 precinct locations to 2010 precincts
    with sqlite3.connect(vpap.va_manualprecinctlocations_filename(2017)) as db:
        for id, name, address, lat, lon in db.execute('SELECT * FROM manual_precinct_locations'):
            precinct_map[id] = index.region_for_point((float(lon), float(lat)))
            #centroids[id] = (float(lon), float(lat))
            assert precinct_map[id] is not None
    
    return precinct_map

def load_vtd_populations():
    vtd_pop_dict = {}
    with sqlite3.connect(census.sf1_geo_state_db_filename('Virginia', 'va', '51')) as db:
        for state, county, vtd, pop100 in db.execute(
            'SELECT STATE, COUNTY, VTD, SUM(POP100) AS POP100 FROM sf1 WHERE SUMLEV = "101" GROUP BY STATE, COUNTY, VTD'
        ):
            vtd_pop_dict[(state, county, int(vtd))] = float(pop100)
    return vtd_pop_dict
