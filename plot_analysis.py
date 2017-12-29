#!/usr/bin/env python

import os
import sys
import shapefile
import matplotlib
from matplotlib import pyplot
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
from matplotlib.figure import Figure
from matplotlib.colorbar import ColorbarBase
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy
import sqlite3

SCRIPT_DIR = os.path.dirname(__file__)
sys.path.append(SCRIPT_DIR)

import census
import vpap
import va2017results
from geometry import *
from geography import *
from districting import *

N_DISTRICTS = 100
ASPECT_RATIO = spherical_distance((-79, 37), (-79, 38)) / spherical_distance((-79, 37), (-78, 37))
RESULTS_DIR = os.path.join(SCRIPT_DIR, 'results')


def main():
    check_files()
    
    precincts = list(iter_precincts_2010())
    precinct_attrs = load_precinct_attributes()
    
    for subdir in os.listdir(RESULTS_DIR):
        results_dir = os.path.join(RESULTS_DIR, subdir)
        db_filename = os.path.join(results_dir, 'results.sqlite')
        if os.path.exists(db_filename):
            D = load_districting(precincts, db_filename)
            district_attrs = load_district_attributes(db_filename)
            
            plot_precincts_popdens(D, results_dir)
            plot_precincts_dempower(D, precinct_attrs, results_dir)
            plot_districts_dempower(D, district_attrs, results_dir)
            plot_districts_winner(D, district_attrs, results_dir)
            plot_districts_distance(D, district_attrs, results_dir)
            plot_districts_pop(D, results_dir)
            

def load_districting(precincts, db_filename):
    region_district_map = {}
    with sqlite3.connect(db_filename) as db:
        for district_id, precinct_id in db.execute('SELECT * FROM district_precincts'):
            region_district_map[precincts[precinct_id]] = district_id
    return Districting(precincts, region_district_map = region_district_map)

def load_precinct_attributes():
    precinct_attrs = {}
    with sqlite3.connect(os.path.join(SCRIPT_DIR, 'data', 'va2017results', 'va2017results.sqlite')) as db:
        for precinct_id, population, dem_votes, rep_votes in db.execute('SELECT * FROM votes_2010precincts'):
            precinct_attrs[precinct_id] = {
                'population': population,
                'dem_votes': dem_votes,
                'rep_votes': rep_votes,
                'dem_power': 0.5 if dem_votes == 0 and rep_votes == 0 else dem_votes / float(dem_votes + rep_votes)
            }
    return precinct_attrs

def load_district_attributes(db_filename):
    district_attrs = {}
    with sqlite3.connect(db_filename) as db:
        for id, population, avg_distance, dem_votes, rep_votes in db.execute('SELECT id, population, avg_distance, dem_votes, rep_votes FROM districts'):
            district_attrs[id] = {
                'population': population,
                'avg_distance' : avg_distance,
                'dem_votes': dem_votes,
                'rep_votes': rep_votes,
                'dem_power': 0.5 if dem_votes == 0 and rep_votes == 0 else dem_votes / float(dem_votes + rep_votes)
            }
    return district_attrs

def plot_districts_dempower(D, district_attrs, results_dir):
    fig, ax = begin_figure()
    precincts = D.regions
    cmap_name = 'bwr'
    valrange = (0.0, 1.0)
    name = 'districts_dempower'
    reppower = {
        d : 1.0 - district_attrs[d]['dem_power'] for d in district_attrs
    }
    plot_district_polygons(fig, ax, D, reppower, cmap_name, valrange)
    end_figure(fig, ax, results_dir, name)
    plot_colorbar(cmap_name, valrange, results_dir, name)

def plot_districts_winner(D, district_attrs, results_dir):
    fig, ax = begin_figure()
    precincts = D.regions
    name = 'districts_winner'
    cmap_name = 'bwr'
    valrange = (0.0, 1.0)
    winner = {
        d : 0.1 if district_attrs[d]['dem_power'] > 0.5 else 0.9
        for d in district_attrs
    }
    plot_district_polygons(fig, ax, D, winner, cmap_name, valrange)
    plot_district_boundaries(fig, ax, D)
    end_figure(fig, ax, results_dir, name)
    plot_colorbar(cmap_name, valrange, results_dir, name)

def plot_districts_distance(D, district_attrs, results_dir):
    fig, ax = begin_figure()
    precincts = D.regions
    name = 'districts_distance'
    cmap_name = 'Oranges'
    valrange = (0, 100000)
    avg_dist = {
        d : attrs['avg_distance']
        for d, attrs in district_attrs.items()
    }
    plot_district_polygons(fig, ax, D, avg_dist, cmap_name, valrange)
    end_figure(fig, ax, results_dir, name)
    plot_colorbar(cmap_name, valrange, results_dir, name)

def plot_precincts_dempower(D, precinct_attrs, results_dir):
    fig, ax = begin_figure()
    precincts = D.regions
    name = 'precincts_dempower'
    cmap_name = 'bwr'
    valrange = (0.0, 1.0)
    reppower = [
        0.5 if not precinct.id in precinct_attrs else 1.0 - precinct_attrs[precinct.id]['dem_power']
        for precinct in precincts
    ]
    plot_precinct_polygons(fig, ax, precincts, reppower, cmap_name, valrange)
    plot_precinct_boundaries(fig, ax, precincts)
    plot_district_boundaries(fig, ax, D)
    end_figure(fig, ax, results_dir, name)
    plot_colorbar(cmap_name, valrange, results_dir, name)

def plot_districts_pop(D, results_dir):
    fig, ax = begin_figure()
    name = 'districts_pop'
    cmap_name = 'Greens'
    valrange = (35000, 145000)
    name = 'districts_pop'
    pop = {d : D.district_population(d) for d in D.districts}
    plot_district_polygons(fig, ax, D, pop, cmap_name, valrange)
    end_figure(fig, ax, results_dir, name)
    plot_colorbar(cmap_name, valrange, results_dir, name)

def plot_precincts_popdens(D, results_dir):
    precincts = D.regions
    fig, ax = begin_figure()
    name = 'precincts_popdens'
    cmap_name = 'Greens'
    logpopdens = [
        log((1.0 + precinct.population) / numpy.abs(precinct.area))
        for precinct in precincts
    ]
    valrange = (min(logpopdens), max(logpopdens))
    plot_precinct_polygons(fig, ax, precincts, logpopdens, cmap_name, valrange)
    plot_precinct_boundaries(fig, ax, precincts)
    plot_district_boundaries(fig, ax, D)
    end_figure(fig, ax, results_dir, name)
    plot_colorbar(cmap_name, valrange, results_dir, name)

def get_colors(x, cmap, valrange):
    xmin = min(x) if valrange is None else valrange[0]
    xmax = max(x) if valrange is None else valrange[1]
    return cmap([(xi - xmin) / (xmax - xmin) for xi in x])

def begin_figure():
    fig = pyplot.figure(figsize = (16, 8), dpi = 600)
    ax = pyplot.gca()
    
    return fig, ax

def end_figure(fig, ax, results_dir, name):
    ax.axes.set_aspect(ASPECT_RATIO)
    ax.autoscale()
    pyplot.axis('off')
    pyplot.savefig(os.path.join(results_dir, f'{name}.png'), bbox_inches='tight')
    pyplot.close(fig)

def plot_colorbar(cmap_name, valrange, results_dir, name):
    fig = pyplot.figure(figsize = (2, 8), dpi = 300)
    ax = pyplot.axes()
    
    colorbar = ColorbarBase(
        ax,
        norm = Normalize(vmin = valrange[0], vmax = valrange[1]),
        cmap = pyplot.get_cmap(cmap_name)
    )
    
    pyplot.subplots_adjust(left = 0.25, right = 0.5)
    ax.autoscale()
    pyplot.savefig(os.path.join(results_dir, f'{name}_colorbar.png'), bbox_inches='tight')
    
    pyplot.close(fig)

def plot_district_boundaries(fig, ax, D):
    lines = []
    for d in D.districts:
        for segment in D.district_polygon(d):
            lines.append(segment)
    
    lc = LineCollection(lines, linewidth = 1.0, color = (0.1, 0.1, 0.1))
    ax.add_collection(lc)

def plot_district_polygons(fig, ax, D, values, cmap_name, valrange, **kwargs):
    values = [
        values[D.region_district(precinct)]
        for precinct in D.regions
    ]
    plot_precinct_polygons(fig, ax, D.regions, values, cmap_name, valrange, **kwargs)

def plot_precinct_polygons(fig, ax, precincts, values, cmap_name, valrange, **kwargs):
    cmap = pyplot.get_cmap(cmap_name)
    colors = get_colors(values, cmap, valrange)
    
    patches = []
    for precinct in precincts:
        polygon = Polygon(precinct.polygon[:-1], closed = True)
        patches.append(polygon)
    collection = PatchCollection(patches, facecolors = colors, edgecolors = colors, **kwargs)
#     collection.set_color(colors)
    ax.add_collection(collection)

def plot_precinct_boundaries(fig, ax, precincts, **kwargs):
    patches = []
    for precinct in precincts:
        polygon = Polygon(precinct.polygon[:-1], closed = True)
        patches.append(polygon)
    collection = PatchCollection(patches, facecolors = 'None', edgecolor = 'darkgray', linewidth = 0.2, **kwargs)
    ax.add_collection(collection)

def check_files():
    # Make sure VA 2017 results database is ready
    va2017results.results_db_filename()

if __name__ == '__main__':
    main()
