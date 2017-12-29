#!/usr/bin/env python

import os
import sys
import shapefile

SCRIPT_DIR = os.path.dirname(__file__)
sys.path.append(SCRIPT_DIR)

import census
import vpap
import va2017results
from geometry import *
from districting import *

N_DISTRICTS = 100

RESULTS_DIR = os.path.join(SCRIPT_DIR, 'results')

def main():
    check_files()
    
    precincts2010 = list(iter_precincts_2010())
    precinct_map = load_precinct_map(precincts2010)
    votes_by_precinct = load_votes_by_precinct(precinct_map)
    
    D0 = Districting(precincts2010)
    
    # Greedy agglomeration step: start with each precinct in its own district,
    # and join neighbors according to some criterion until we have N_DISTRICTS districts.
    for agglom_func, swap_func in (
        # Minimize total distance between people within district.
        (agglomerate_mindist, swap_mindist),
        
        # Maximize population balance entropy
        (agglomerate_balance, swap_balance),
        
        # Stepwise minimize average rank of min. dist, max. entropy
        (agglomerate_hybrid, swap_hybrid),
    ):
        print('Agglomerating using {}...'.format(agglom_func.__name__))
        D_agglomerated = agglom_func(D0, N_DISTRICTS)
        D_agglomerated.verify()
        write_districting_results(D_agglomerated, agglom_func.__name__, votes_by_precinct)
        
        print('Swapping using {}...'.format(swap_func.__name__))
        D_swapped = swap_func(D_agglomerated)
        D_swapped.verify()
        write_districting_results(D_swapped, '{}+{}'.format(
            agglom_func.__name__,
            swap_func.__name__
        ), votes_by_precinct)

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

def load_votes_by_precinct(precinct_map):
    # Count votes in 2010 precincts
    votes_by_precinct = {}
    for precinct in precinct_map.values():
        votes_by_precinct[precinct] = {}
    with sqlite3.connect(va2017results.results_db_filename()) as db:
        for precinct_id, party, votes in db.execute('''
            SELECT precinct_id, party, votes FROM votes, candidates WHERE candidate_id = candidates.id
        '''):
            precinct = precinct_map[precinct_id]
            if not party in votes_by_precinct[precinct]:
                votes_by_precinct[precinct][party] = 0
            votes_by_precinct[precinct][party] += int(votes)
    return votes_by_precinct

def map_votes_to_districts(votes_by_precinct, D):
    votes_by_district = {}
    for d in D.districts:
        votes_by_district[d] = {}
    
    for precinct, votes_by_party in votes_by_precinct.items():
        district = D.region_district(precinct)
        for party, votes in votes_by_party.items():
            if not party in votes_by_district[district]:
                votes_by_district[district][party] = 0
            votes_by_district[district][party] += int(votes)
    
    return votes_by_district

def iter_precinct_shaperecords():
    sf = shapefile.Reader(census.vtd_state_shapefiles_filename('Virginia', 'va', '51'))
    for sr in sf.shapeRecords():
        yield sr.shape, {field[0] : value for field, value in zip(sf.fields[1:], sr.record)}
    sf.close()

def load_vtd_populations():
    vtd_pop_dict = {}
    with sqlite3.connect(census.sf1_geo_state_db_filename('Virginia', 'va', '51')) as db:
        for state, county, vtd, pop100 in db.execute(
            'SELECT STATE, COUNTY, VTD, SUM(POP100) AS POP100 FROM sf1 WHERE SUMLEV = "101" GROUP BY STATE, COUNTY, VTD'
        ):
            vtd_pop_dict[(state, county, int(vtd))] = float(pop100)
    return vtd_pop_dict

def write_districting_results(D, name, votes_by_precinct):
    votes_by_district = map_votes_to_districts(votes_by_precinct, D)
    
    results_subdir = os.path.join(RESULTS_DIR, name)
    try:
        os.makedirs(results_subdir)
    except:
        pass
    
    db_filename = os.path.join(results_subdir, 'results.sqlite')
    with sqlite3.connect(db_filename) as db:
        c = db.cursor()
        
        c.execute('CREATE TABLE district_precincts (district_id, precinct_id)')
        c.execute('CREATE INDEX idx_district_precincts ON district_precincts (district_id, precinct_id)')
        
        for d in D.districts:
            for p in D.district_regions(d):
                c.execute('INSERT INTO district_precincts VALUES (?,?)', (d, p.id))
        
        c.execute('CREATE TABLE districts (id, population, area, avg_distance, dem_votes, rep_votes)')
        c.execute('CREATE INDEX idx_districts ON districts (id, dem_votes, rep_votes)')
        
        for d, votes_by_party in votes_by_district.items():
            dem_votes = votes_by_party['Democratic'] if 'Democratic' in votes_by_party else 0
            rep_votes = votes_by_party['Republican'] if 'Republican' in votes_by_party else 0
            c.execute(
                'INSERT INTO districts VALUES (?,?,?,?,?,?)',
                (d, D.district_population(d), D.district_area(d), D.district_average_distance(d), dem_votes, rep_votes)
            )
        
        c.execute('CREATE TABLE power (dem_seats, rep_seats, dem_pop, rep_pop, dem_power, rep_power)')
        
        dem_seats, = c.execute('SELECT COUNT(*) FROM districts WHERE dem_votes > rep_votes').fetchone()
        rep_seats, = c.execute('SELECT COUNT(*) FROM districts WHERE rep_votes > dem_votes').fetchone()
        
        dem_pop, = c.execute('''
            SELECT SUM(population) FROM districts
            WHERE dem_votes > rep_votes
        ''').fetchone()
        rep_pop, = c.execute('''
            SELECT SUM(population) FROM districts
            WHERE rep_votes > dem_votes
        ''').fetchone()
        
        dem_power = dem_pop / float(dem_pop + rep_pop)
        rep_power = rep_pop / float(dem_pop + rep_pop)
        
        c.execute(
            'INSERT INTO power VALUES (?,?,?,?,?,?)',
            (dem_seats, rep_seats, dem_pop, rep_pop, dem_power, rep_power)
        )

def check_files():
    # Make sure census database is ready
    census.sf1_geo_state_db_filename('Virginia', 'va', '51')
    
    # Make sure VA precinct shapefiles from 2010 and 2016 are ready
    census.vtd_state_shapefiles_filename('Virginia', 'va', '51')
    vpap.va_precinctmaps_filename(2016)
    vpap.va_manualprecinctlocations_filename(2017)
    
    # Make sure VA 2017 results database is ready
    va2017results.results_db_filename()

if __name__ == '__main__':
    main()
