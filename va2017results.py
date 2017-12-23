import os
from collections import OrderedDict
import json
import sqlite3
from download import download_url

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data', 'va2017results')
URL_PREFIX = 'http://results.elections.virginia.gov/vaelections/2017 November General/Json'

LOCALITIES = [
    'ACCOMACK COUNTY',
    'ALBEMARLE COUNTY',
    'ALEXANDRIA CITY',
    'ALLEGHANY COUNTY',
    'AMELIA COUNTY',
    'AMHERST COUNTY',
    'APPOMATTOX COUNTY',
    'ARLINGTON COUNTY',
    'AUGUSTA COUNTY',
    'BATH COUNTY',
    'BEDFORD COUNTY',
    'BLAND COUNTY',
    'BOTETOURT COUNTY',
    'BRISTOL CITY',
    'BRUNSWICK COUNTY',
    'BUCHANAN COUNTY',
    'BUCKINGHAM COUNTY',
    'BUENA VISTA CITY',
    'CAMPBELL COUNTY',
    'CAROLINE COUNTY',
    'CARROLL COUNTY',
    'CHARLES CITY COUNTY',
    'CHARLOTTE COUNTY',
    'CHARLOTTESVILLE CITY',
    'CHESAPEAKE CITY',
    'CHESTERFIELD COUNTY',
    'CLARKE COUNTY',
    'COLONIAL HEIGHTS CITY',
    'COVINGTON CITY',
    'CRAIG COUNTY',
    'CULPEPER COUNTY',
    'CUMBERLAND COUNTY',
    'DANVILLE CITY',
    'DICKENSON COUNTY',
    'DINWIDDIE COUNTY',
    'EMPORIA CITY',
    'ESSEX COUNTY',
    'FAIRFAX CITY',
    'FAIRFAX COUNTY',
    'FALLS CHURCH CITY',
    'FAUQUIER COUNTY',
    'FLOYD COUNTY',
    'FLUVANNA COUNTY',
    'FRANKLIN CITY',
    'FRANKLIN COUNTY',
    'FREDERICK COUNTY',
    'FREDERICKSBURG CITY',
    'GALAX CITY',
    'GILES COUNTY',
    'GLOUCESTER COUNTY',
    'GOOCHLAND COUNTY',
    'GRAYSON COUNTY',
    'GREENE COUNTY',
    'GREENSVILLE COUNTY',
    'HALIFAX COUNTY',
    'HAMPTON CITY',
    'HANOVER COUNTY',
    'HARRISONBURG CITY',
    'HENRICO COUNTY',
    'HENRY COUNTY',
    'HIGHLAND COUNTY',
    'HOPEWELL CITY',
    'ISLE OF WIGHT COUNTY',
    'JAMES CITY COUNTY',
    'KING & QUEEN COUNTY',
    'KING GEORGE COUNTY',
    'KING WILLIAM COUNTY',
    'LANCASTER COUNTY',
    'LEE COUNTY',
    'LEXINGTON CITY',
    'LOUDOUN COUNTY',
    'LOUISA COUNTY',
    'LUNENBURG COUNTY',
    'LYNCHBURG CITY',
    'MADISON COUNTY',
    'MANASSAS CITY',
    'MANASSAS PARK CITY',
    'MARTINSVILLE CITY',
    'MATHEWS COUNTY',
    'MECKLENBURG COUNTY',
    'MIDDLESEX COUNTY',
    'MONTGOMERY COUNTY',
    'NELSON COUNTY',
    'NEW KENT COUNTY',
    'NEWPORT NEWS CITY',
    'NORFOLK CITY',
    'NORTHAMPTON COUNTY',
    'NORTHUMBERLAND COUNTY',
    'NORTON CITY',
    'NOTTOWAY COUNTY',
    'ORANGE COUNTY',
    'PAGE COUNTY',
    'PATRICK COUNTY',
    'PETERSBURG CITY',
    'PITTSYLVANIA COUNTY',
    'POQUOSON CITY',
    'PORTSMOUTH CITY',
    'POWHATAN COUNTY',
    'PRINCE EDWARD COUNTY',
    'PRINCE GEORGE COUNTY',
    'PRINCE WILLIAM COUNTY',
    'PULASKI COUNTY',
    'RADFORD CITY',
    'RAPPAHANNOCK COUNTY',
    'RICHMOND CITY',
    'RICHMOND COUNTY',
    'ROANOKE CITY',
    'ROANOKE COUNTY',
    'ROCKBRIDGE COUNTY',
    'ROCKINGHAM COUNTY',
    'RUSSELL COUNTY',
    'SALEM CITY',
    'SCOTT COUNTY',
    'SHENANDOAH COUNTY',
    'SMYTH COUNTY',
    'SOUTHAMPTON COUNTY',
    'SPOTSYLVANIA COUNTY',
    'STAFFORD COUNTY',
    'STAUNTON CITY',
    'SUFFOLK CITY',
    'SURRY COUNTY',
    'SUSSEX COUNTY',
    'TAZEWELL COUNTY',
    'VIRGINIA BEACH CITY',
    'WARREN COUNTY',
    'WASHINGTON COUNTY',
    'WAYNESBORO CITY',
    'WESTMORELAND COUNTY',
    'WILLIAMSBURG CITY',
    'WINCHESTER CITY',
    'WISE COUNTY',
    'WYTHE COUNTY',
    'YORK COUNTY',
]

def results_db_filename():
    filename = os.path.join(DATA_DIR, 'va2017results.sqlite')
    
    if not os.path.exists(filename):
        
        with sqlite3.connect(filename) as db:
            c = db.cursor()
        
            c.execute('''CREATE TABLE candidates (
                id INTEGER PRIMARY KEY, seat_id, name, party
            )''')
            c.execute('''CREATE TABLE localities (
                id PRIMARY KEY, name
            )''')
            c.execute('''CREATE TABLE precincts (
                id PRIMARY KEY, name, locality_id
            )''')
            c.execute('CREATE INDEX idx_candidates_seat_id_name ON candidates (seat_id, name)')
            c.execute('CREATE INDEX idx_candidates_party ON candidates (seat_id, party)')
            c.execute('''CREATE TABLE votes (
                precinct_id, seat_id, candidate_id, votes
            )''')
            c.execute('CREATE INDEX idx_votes_precinct_id ON votes (precinct_id)')
            
            for locality in LOCALITIES:
                locality_data = locality_aggregate_data(locality)
                for race_obj in locality_data['Races']:
                    contest = race_obj['RaceName']
                    if contest.startswith('Member House of Delegates'):
                        load_contest_into_db(locality, contest, c)
    
    assert os.path.exists(filename)
    return filename

def load_contest_into_db(locality, contest, c):
    data = locality_contest_byprecinct_data(locality, contest)
    
    locality_id = data['Locality']['LocalityCode']
    assert data['RaceName'].startswith('Member House of Delegates (')
    assert len(data['RaceName']) == len('Member House of Delegates (XXX)')
    seat_id = data['RaceName'][-4:-1]
    
    c.execute(
        'INSERT OR IGNORE INTO localities (id, name) VALUES (?, ?)',
        (locality_id, data['Locality']['LocalityName'])
    )
    
    print('Gathering data for locality {}, seat {}'.format(locality_id, seat_id))
    for precinct in data['Precincts']:
        if precinct['PrecinctName'].startswith('#'):
            print('Skipping unidentifiable precinct {}'.format(precinct['PrecinctName']))
        else:
            precinct_id = '51{}{}'.format(
                locality_id, precinct['PrecinctName'][:3]
            )
            try:
                precinct_name = precinct['PrecinctName'][:-6].split('-')[1].strip()
            except Exception as e:
                precinct_name = precinct['PrecinctName'][4:-6]
            c.execute(
                'INSERT OR IGNORE INTO precincts (id, name, locality_id) VALUES (?,?,?)',
                (precinct_id, precinct_name, locality_id)
            )
            
            print('Processing precinct {}, {}'.format(precinct_id, precinct_name))
            for candidate in precinct['Candidates']:
                name = candidate['BallotName']
                get_candidate_query = 'SELECT candidate_id FROM candidates WHERE seat_id = ? AND name = ?'

                try:
                    candidate_id, = c.execute(get_candidate_query, [seat_id, name]).fetchone()
                except:
                    c.execute(
                        'INSERT INTO candidates (seat_id, name, party) VALUES (?,?,?)',
                        [seat_id, name, candidate['PoliticalParty']]
                    )
                    candidate_id = c.lastrowid
                    print(candidate_id)
                c.execute(
                    'INSERT INTO votes (precinct_id, seat_id, candidate_id, votes) VALUES (?, ?, ?, ?)',
                    [precinct_id, seat_id, candidate_id, candidate['Votes']]
                )

def general_assembly_aggregate_filename():
    filename = os.path.join(DATA_DIR, 'general_assembly_aggregate.json')
    if not os.path.exists(filename):
        download_url(
            '{}/GeneralAssembly.json'.format(URL_PREFIX),
            filename
        )
    return filename

def locality_aggregate_data(locality):
    with open(locality_aggregate_filename(locality)) as f:
        return json.load(f)

def locality_aggregate_filename(locality):
    filename = os.path.join(DATA_DIR, 'locality_aggregate', locality, 'Index.json')
    if not os.path.exists(filename):
        download_url(
            '{}/Locality/{locality}/Index.json'.format(URL_PREFIX),
            filename
        )
    return filename

def locality_contest_byprecinct_data(locality, contest):
    with open(locality_contest_byprecinct_filename(locality, contest)) as f:
        return json.load(f)

def locality_contest_byprecinct_filename(locality, contest):
    filename = os.path.join(DATA_DIR, 'locality_byprecinct', locality, '{}.json'.format(contest))
    if not os.path.exists(filename):
        download_url(
            '{}/Locality/{}/{}.json'.format(URL_PREFIX, locality, contest),
            filename
        )
    return filename
