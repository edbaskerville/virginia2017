import os
import sys
SCRIPT_DIR = os.path.dirname(__file__)
sys.path.append(SCRIPT_DIR)
from geometry import *
import sqlite3
from collections import deque
from sortedcontainers import SortedListWithKey

class Region:
    def __init__(self, id, polygon, population = 0, record = None):
        self.id = id
        self.polygon = polygon
        self.population = population
        self.record = record
        
        try:
            assert self.polygon[0] == self.polygon[-1]
        except Exception as e:
            print(self.polygon[0], self.polygon[-1])
            raise e
        
        self._centroid = None
        self._area = None
        
        self.xmin = min(x for x, y in polygon)
        self.xmax = max(x for x, y in polygon)
        self.ymin = min(y for x, y in polygon)
        self.ymax = max(y for x, y in polygon)
    
    def iter_segments(self):
        V = self.polygon
        for i in range(len(V) - 1):
            yield (V[i], V[i+1])
    
    def compute_total_distance_per_capita(self, regions):
        return sum(
            region.population * spherical_distance(d2r(self.centroid), d2r(region.centroid))
            for region in regions if region is not self
        )
    
    @property
    def centroid(self):
        if self._centroid is None:
            self._centroid = polygon_centroid(self.polygon, self.area)
        return self._centroid
    
    @property
    def area(self):
        if self._area is None:
            self._area = polygon_area(self.polygon)
        return self._area
    
    def contains_point(self, p):
        x, y = p
        
        # Fast bounds check
        if x < self.xmin:
            return False
        if x > self.xmax:
            return False
        if y < self.ymin:
            return False
        if y > self.ymax:
            return False
        
        return polygon_contains_point(self.polygon, p)
    
    def __lt__(self, other):
        return self.id < other.id
    
    def __eq__(self, other):
        if isinstance(other, Region):
            return self.id == other.id
        return False
    
    def __hash__(self):
        return int(self.id)
    
    def __repr__(self):
        return 'Region({})'.format(self.id)

class RegionIndex:
    def __init__(self, regions):
        self.regions = regions
        self.by_lon = SortedListWithKey(regions, key = lambda r : r.centroid[0])
        self.by_lat = SortedListWithKey(regions, key = lambda r : r.centroid[1])
    
    def get_closest(self, l, x):
        index = l.bisect_key(x)
        if index == 0:
            return [l[0]]
        if index == len(l):
            return [l[-1]]
        return [l[index - 1], l[index]]
    
    # TODO: this could probably theoretically fail under certain circumstances
    def region_for_point(self, point):
        lon, lat = point
        
        # Identify distance to closest region
        d = min(
            euclidean_distance(r.centroid, point)
            for r in (
                self.get_closest(self.by_lon, lon) +
                self.get_closest(self.by_lat, lat)
            )
        )
        
        # Go through all regions whose centroid is in the 2*d X 2*d box around this point
        #lon_delta = sin(lat * pi / 180.0) * d * 180.0 / pi
        test_regions = set(
            self.by_lon.irange_key(lon - d, lon + d, inclusive = (True, True))
        ) & set(
            self.by_lat.irange_key(lat - d, lat + d, inclusive = (True, True))
        )
        
        for r in test_regions:
            if r.contains_point(point):
                return r
        
        # If we don't find one, return the one with the closest centroid
        assert len(test_regions) > 0
        return min(test_regions, key = lambda r : euclidean_distance(r.centroid, point))

def sorted_tuple(*objs):
    return tuple(sorted(objs))

class Graph:
    def __init__(self, vertices, edges):
        self.vertices = set(vertices)
        self.edges = set(sorted_tuple(*edge) for edge in edges)
        
        self.neighbor_map = {v: set() for v in self.vertices}
        for v, w in self.edges:
            self.neighbor_map[v].add(w)
            self.neighbor_map[w].add(v)
        
        for v, neighbors in self.neighbor_map.items():
            assert len(neighbors) > 0
    
    @property
    def vertex_count(self):
        return len(self.vertices)
    
    @property
    def edge_count(self):
        return len(self.edges)
    
    def neighbors(self, v):
        return self.neighbor_map[v]
    
    def contains_vertices(self, v):
        return v in self.nodes
    
    def contains_edge(self, v, w):
        return sorted_tuple(v, w) in self.edges

def make_region_graph(regions):
    seg_reg_map = {}
    for i, r in enumerate(regions):
        for segment in r.iter_segments():
            segment = sorted_tuple(*segment)
            if not segment in seg_reg_map:
                seg_reg_map[segment] = set()
            seg_reg_map[segment].add(r)
    
    def iterate_edges():
        for segment, segment_regions in seg_reg_map.items():
            segment_regions = list(segment_regions)
            if len(segment_regions) == 2:
                r1, r2 = segment_regions
                yield r1, r2
            else:
                for i in range(len(segment_regions) - 1):
                    for j in range(i+1, len(segment_regions)):
                        yield segment_regions[i], segment_regions[j]
    
    return Graph(regions, iterate_edges())

class Districting:
    def __init__(
        self, regions, region_district_map = None,
        _population = None,
        _districts = None,
        _region_graph = None, 
        _district_region_map = None,
        _boundary_regions = None,
        _neighbor_districts = None,
        _total_distance = None,
        _merge_distance_delta = None,
        _move_distance_deltas = None
    ):
        self.regions = list(regions)
        self._population = _population
        
        if region_district_map is None:
            self.region_district_map = {r : r.id for r in regions}
        else:
            self.region_district_map = region_district_map
        
        self._districts = _districts
        
        self._region_graph = _region_graph
        self._district_region_map = _district_region_map
        
        self._boundary_regions = _boundary_regions
        self._neighbor_districts = _neighbor_districts
        
        self._total_distance = _total_distance
        self._merge_distance_delta = _merge_distance_delta
        
        self._move_distance_deltas = _move_distance_deltas
        
        self._district_populations = None
        
        self._set = None
        self._colors = None
    
    def verify(self):
        self.verify_contiguity()
    
    def verify_contiguity(self):
        for d in self.districts:
            # Make sure we hit every region in the district
            visited = set()
            q = deque()
            q.append(next(iter(self.district_regions(d))))
            while len(q) > 0:
                r = q.popleft()
                visited.add(r)
                for rn in self.neighbor_regions(r):
                    if not rn in visited and self.region_district(rn) == d:
                        q.append(rn)
            print(len(self.district_regions(d)), len(visited))
            assert len(self.district_regions(d) & visited) == len(self.district_regions(d))
            assert len(self.district_regions(d) - visited) == 0
    
    @property
    def total_distance(self):
        if self._total_distance is None:
            self._total_distance = sum(
                self.compute_total_distance(d, d) / 2.0 for d in self.districts
            )
        return self._total_distance
    
    def district_average_distance(self, d):
        return self.compute_total_distance(d, d) / (self.district_population(d)**2.0)
    
    def district_area(self, d):
        return sum(abs(r.area) for r in self.district_regions(d))
    
    def merge_distance_delta(self, d1, d2):
        k = sorted_tuple(d1, d2)
        
        if self._merge_distance_delta is None:
            self._merge_distance_delta = {}
        if k not in self._merge_distance_delta:
            self._merge_distance_delta[k] = self.compute_total_distance(d1, d2)
        
        return self._merge_distance_delta[k]
    
    def merge_entropy_delta(self, d1, d2):
        pop1 = self.district_population(d1)
        pop2 = self.district_population(d2)
        
        p1 = pop1 / self.population
        p2 = pop2 / self.population
        
        p = (pop1 + pop2) / self.population
        
        delta = (
            (
                (0.0 if p1 == 0.0 else p1 * log(p1)) +
                (0.0 if p2 == 0.0 else p2 * log(p2))
            ) -
            (
                0.0 if p == 0.0 else p * log(p)
            )
        )
        #print(pop1, pop2, p1, p2, p, delta)
        return delta
    
    def merge_population_below_ceiling(self, d1, d2, k):
        target_pop = self.population / k
        pop1 = self.district_population(d1)
        pop2 = self.district_population(d2)
        return pop1 + pop2 <= target_pop
    
    def compute_total_distance(self, d1, d2):
        return sum(
            r.population * r.compute_total_distance_per_capita(self.district_region_map[d2])
            for r in self.district_region_map[d1]
        )
    
    def merge_merge_distance_delta(self, d1, d2):
        merge_distance_delta = dict(self._merge_distance_delta)
        for dn in self._neighbor_districts[d2]:
            k = sorted_tuple(d2, dn)
            if k in merge_distance_delta:
                del merge_distance_delta[k]
        for dn in self._neighbor_districts[d1] | self._neighbor_districts[d2]:
            if dn != d1 and dn != d2:
                merge_distance_delta[sorted_tuple(d1, dn)] = \
                    self.merge_distance_delta(d1, dn) + self.merge_distance_delta(d2, dn)
        return merge_distance_delta
    
    def merge_districts(self, d1, d2):
        region_district_map = dict(self.region_district_map)
        for r in self.district_region_map[d2]:
            region_district_map[r] = d1
        
        district_region_map = dict(self.district_region_map)
        del district_region_map[d2]
        district_region_map[d1] = self.district_region_map[d1] | self.district_region_map[d2]
        
        neighbor_districts = dict(self._neighbor_districts)
        neighbor_districts[d1] = (self._neighbor_districts[d1] | self._neighbor_districts[d2]) - frozenset((d1, d2))
        assert not d2 in neighbor_districts[d1]
        for dn in self._neighbor_districts[d2]:
            if dn != d1:
                neighbor_districts[dn] = (self._neighbor_districts[dn] | frozenset((d1,))) - frozenset((d2,))
        del neighbor_districts[d2]

        for d, ds in neighbor_districts.items():
            assert d != d2
            if d2 in ds:
                assert d in self._neighbor_districts[d2]
            assert not d2 in ds
        
        districts = [d for d in self._districts if d != d2]
        total_distance = self.total_distance + self.merge_distance_delta(d1, d2)
        
        merge_distance_delta = self.merge_merge_distance_delta(d1, d2)
        
        return Districting(
            self.regions,
            region_district_map = region_district_map,
            _region_graph = self._region_graph,
            _district_region_map = district_region_map,
            _neighbor_districts = neighbor_districts,
            _districts = districts,
            _total_distance = total_distance,
            _population = self._population,
            _merge_distance_delta = merge_distance_delta
        )
    
    def move_region(self, r, dnew):
        return Districting(
            self.regions,
            region_district_map = self.region_district_map_after_move(r, dnew),
            _region_graph = self._region_graph,
            _district_region_map  = self.district_region_map_after_move(r, dnew),
            _districts = self._districts,
            _total_distance = self.total_distance_after_move(r, dnew),
            _boundary_regions = self.boundary_regions_after_move(r, dnew),
            _move_distance_deltas = self.move_distance_deltas_after_move(r, dnew)
        )
    
    def region_district_map_after_move(self, r, dnew):
        rdm = dict(self.region_district_map)
        rdm[r] = dnew
        return rdm
    
    def district_region_map_after_move(self, r, dnew):
        dold = self.region_district_map[r]
        drm = dict(self._district_region_map)
        rset = frozenset((r,))
        drm[dold] = drm[dold] - rset
        drm[dnew] = drm[dnew] | rset
        return drm
    
    def total_distance_after_move(self, r, dnew):
        return self.total_distance + self.move_distance_delta(r, dnew)
    
    def move_distance_delta(self, r, dnew):
        dold = self.region_district(r)
        assert dnew != dold
        
        if self._move_distance_deltas is None:
            self._move_distance_deltas = {}
        if self._move_distance_deltas is None:
            self._move_distance_deltas = {}
        if not r in self._move_distance_deltas:
            self._move_distance_deltas[r] = {}
        
        if not dnew in self._move_distance_deltas[r]:
            delta_old = -r.population * r.compute_total_distance_per_capita(self.district_regions(dold))
            delta_new = r.population * r.compute_total_distance_per_capita(self.district_regions(dnew))
            self._move_distance_deltas[r][dnew] = delta_old + delta_new
        
        return self._move_distance_deltas[r][dnew]
    
    @property
    def population_entropy(self):
        p = [self.district_population(d) / self.population for d in self.districts]
        return -sum(
            (0.0 if p_i == 0.0 else p_i * log(p_i))
            for p_i in p
        )
    
    def move_entropy_delta(self, r, d2):
        d1 = self.region_district(r)
        pop_d1_old = self.district_population(d1)
        pop_d2_old = self.district_population(d2)
        pop_d1_new = pop_d1_old - r.population
        pop_d2_new = pop_d2_old + r.population
        
        p_d1_old = pop_d1_old / self.population
        p_d2_old = pop_d2_old / self.population
        
        p_d1_new = pop_d1_new / self.population
        p_d2_new = pop_d2_new / self.population
        
        return (
            (p_d1_old * log(p_d1_old) + p_d2_old * log(p_d2_old)) -
            (p_d1_new * log(p_d1_new) + p_d2_new * log(p_d2_new))
        )
    
    def move_distance_deltas_after_move(self, r, dnew):
        assert self._move_distance_deltas is not None
        
        dold = self.region_district(r)
        
        # Overzealous invalidation: optimize if not fast enough
        mdd = dict(self._move_distance_deltas)
        for rd in self.district_regions(dold):
            if rd in mdd:
                del mdd[rd]
        for rd in self.district_regions(dnew):
            if rd in mdd:
                del mdd[rd]
        for rn in self.neighbor_regions(r):
            if rn in mdd:
                del mdd[rn]
        
        return mdd
    
    def boundary_regions_after_move(self, r, dnew):
        dold = self.region_district_map[r]
        br = dict(self._boundary_regions)
        
        rset = frozenset((r,))
        if dold in br:
            br[dold] = br[dold] - rset
        if dnew in br:
            br[dnew] = br[dnew] | rset
        
        for rn in self.neighbor_regions(r):
            dn = self.region_district(rn)
            if dold in br:
                if dn == dold and not rn in br[dold]:
                    br[dold].add(rn)
            if dnew in br:
                if dn == dnew and rn in br[dnew] and not self.is_boundary_region(rn, r_to_ignore = r):
                    br[dnew].remove(rn)
        
        return br
    
    def is_boundary_region(self, r, r_to_ignore = None):
        d = self.region_district(r)
        d_rs = self.district_regions(d)
        for rn in self.neighbor_regions(r):
            if r != r_to_ignore and rn not in d_rs:
                return True
        return False
    
    def region_district(self, r):
        return self.region_district_map[r]
    
    def neighbor_regions(self, r):
        return self.region_graph.neighbors(r)
    
    def iter_all_merges(self):
        for d1 in self.districts:
            for d2 in self.neighbor_districts(d1):
                if d1 < d2:
                    yield (d1, d2)
    
    def iter_all_moves(self):
        for d in self.districts:
            if self.region_count(d) > 1:
                for r in self.boundary_regions(d):
                    for dn in set(self.region_district(rn) for rn in self.neighbor_regions(r)):
                        if dn != d and not self.removal_disconnects_district(r):
                            yield (r, dn)
    
    def removal_disconnects_district(self, region):
        d = self.region_district(region)
        q = deque()
        visited = set()
        for r in self.neighbor_regions(region):
            if self.region_district(r) == d:
                q.append(r)
                break
        while len(q) > 0:
            r = q.popleft()
            visited.add(r)
            for rn in self.neighbor_regions(r):
                if rn != region and not rn in visited and self.region_district(rn) == d:
                    q.append(rn)
        return len(visited) < len(self.district_regions(d)) - 1
    
    def districts_are_neighbors(self, d1, d2):
        return d2 in self.neighbor_districts(d1)
    
    def boundary_regions(self, d):
        if self._boundary_regions is None:
            self._boundary_regions = {}
        if not d in self._boundary_regions:
            self._boundary_regions[d] = set()
            for r in self.district_region_map[d]:
                for rn in self.region_graph.neighbors(r):
                    if rn not in self.district_region_map[d]:
                        self._boundary_regions[d].add(r)
        return self._boundary_regions[d]
    
    def neighbor_districts(self, d):
        if self._neighbor_districts is None:
            self._neighbor_districts = {}
            for di in self.districts:
                self._neighbor_districts[di] = set()
                for r in self.boundary_regions(di):
                    for rn in self.region_graph.neighbors(r):
                        if rn not in self.district_region_map[di]:
                            self._neighbor_districts[di].add(self.region_district_map[rn])
        return self._neighbor_districts[d]
    
    @property
    def colors(self):
        if self._colors is None:
            colors = {}
            for d in self.districts:
                try:
                    colors[d] = max(
                        colors[dn] for dn in self.neighbor_districts(d)
                        if dn in colors
                    ) + 1
                except:
                    colors[d] = 0
            self._colors = colors
        return self._colors
    
    def district_polygon(self, d):
        segment_count = {}
        for r in self.district_regions(d):
            for i in range(len(r.polygon) - 1):
                segment = sorted_tuple(*r.polygon[i:i+2])
                if segment in segment_count:
                    segment_count[segment] += 1
                else:
                    segment_count[segment] = 1
        return [segment for segment, count in segment_count.items() if count == 1]
    
    @property
    def region_graph(self):
        if self._region_graph is None:
            self._region_graph = make_region_graph(self.regions)
        return self._region_graph
    
    @property
    def districts(self):
        if self._districts is None:
            self._districts = sorted(self.district_region_map.keys())
        return self._districts
    
    @property
    def district_count(self):
        return len(self.districts)
    
    def region_count(self, d):
        return len(self.district_regions(d))
    
    def district_regions(self, d):
        return self.district_region_map[d]
    
    @property
    def district_region_map(self):
        if self._district_region_map is None:
            self._district_region_map = {}
            for region, district in self.region_district_map.items():
                if not district in self._district_region_map:
                    self._district_region_map[district] = set()
                self._district_region_map[district].add(region)
        return self._district_region_map
    
    @property
    def population(self):
        if self._population is None:
            self._population = sum(r.population for r in self.regions)
        return self._population
    
    def district_population(self, district):
        return self.district_populations[district]
    
    @property
    def district_populations(self):
        if self._district_populations is None:
            self._district_populations = {
                d : sum(r.population for r in self.district_regions(d))
                for d in self.districts
            }
        return self._district_populations
    
    @property
    def set(self):
        if self._set is None:
            self._set = frozenset(
                frozenset(r.id for r in self.district_regions(d))
                for d in self.districts
            )
        return self._set

def inverse_permutation(l):
    l2 = [-1] * len(l)
    for i, li in enumerate(l):
        l2[li] = i
    return l2

def sorted_rank(l, key):
    return inverse_permutation([
        i for i, li in sorted(enumerate(l), key = lambda pair: key(pair[1]))
    ])

def agglomerate_hybrid(D0, k, mindist_weight = 0.5):
    Di = D0
    
    w = mindist_weight
    while Di.district_count > k:
        merges = list(Di.iter_all_merges())
        rank_mindist = sorted_rank(
            merges, key = lambda m : Di.merge_distance_delta(*m)
        )
        rank_balance = sorted_rank(
            merges, key = lambda m : -Di.merge_entropy_delta(*m)
        )
        d1, d2 = min(
            enumerate(merges),
            key = lambda pair: w * rank_mindist[pair[0]] + (1.0 - w) * rank_balance[pair[0]]
        )[1]
        Di = Di.merge_districts(d1, d2)
    return Di
    

def agglomerate_mindist(D0, k):
    Di = D0
    while Di.district_count > k:
        d1, d2 = min(
            Di.iter_all_merges(),
            key = lambda m : Di.merge_distance_delta(*m)
        )
        Di = Di.merge_districts(d1, d2)
    return Di

def agglomerate_balance(D0, k):
    Di = D0
    while Di.district_count > k:
        d1, d2 = max(
            Di.iter_all_merges(),
            key = lambda m : Di.merge_entropy_delta(*m)
        )
        Di = Di.merge_districts(d1, d2)
        print(Di.population_entropy)
    return Di

def swap_mindist(D0):
    Di = D0
    while True:
        r, d = min(
            Di.iter_all_moves(),
            key = lambda m : Di.move_distance_delta(*m)
        )
        delta = Di.move_distance_delta(r, d)
        if delta < 0:
            Di = Di.move_region(r, d)
        else:
            return Di

def swap_balance(D0):
    Di = D0
    while True:
        r, d = max(
            Di.iter_all_moves(),
            key = lambda m : Di.move_entropy_delta(*m)
        )
        delta = Di.move_entropy_delta(r, d)
        if delta > 0:
            Di = Di.move_region(r, d)
            print(Di.population_entropy)
        else:
            return Di

def swap_hybrid(D0, mindist_weight = 0.5):
    w = mindist_weight
    Di = D0
    
    visited = set()
    
    while True:
        moves = list(Di.iter_all_moves())
        rank_mindist = sorted_rank(
            moves, key = lambda m : Di.move_distance_delta(*m)
        )
        rank_balance = sorted_rank(
            moves, key = lambda m : -Di.move_entropy_delta(*m)
        )
        r, d = min(
            enumerate(moves),
            key = lambda pair: w * rank_mindist[pair[0]] + (1.0 - w) * rank_balance[pair[0]]
        )[1]
        distance_delta = Di.move_distance_delta(r, d)
        entropy_delta = Di.move_entropy_delta(r, d)
        if distance_delta < 0 or entropy_delta > 0:
            print(distance_delta, entropy_delta)
            Di = Di.move_region(r, d)
        else:
            break
        
        if Di.set in visited:
            break
        visited.add(Di.set)
    
    return Di
