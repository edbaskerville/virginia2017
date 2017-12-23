from math import asin, sqrt, sin, cos, pi, tan, atan, exp, tan, log

def polygon_centroid(polygon, area = None):
    area6 = 6.0 * (polygon_area(polygon) if area is None else area)
    cx = sum(
        (polygon[i][0] + polygon[i+1][0]) * (
            polygon[i][0] * polygon[i+1][1] - polygon[i+1][0] * polygon[i][1]
        )
        for i in range(len(polygon) - 1)
    ) / area6
    cy = sum(
        (polygon[i][1] + polygon[i+1][1]) * (
            polygon[i][0] * polygon[i+1][1] - polygon[i+1][0] * polygon[i][1]
        )
        for i in range(len(polygon) - 1)
    ) / area6
    
    try:
        assert cx > min(p[0] for p in polygon)
        assert cx < max(p[0] for p in polygon)
        assert cy > min(p[1] for p in polygon)
        assert cy < max(p[1] for p in polygon)
    except Exception as e:
        print(cx, cy)
        pyplot.plot(polygon)
        raise e
    
    return cx, cy

def polygon_area(polygon):
    return 0.5 * sum(
        polygon[i][0] * polygon[i+1][1] - polygon[i+1][0] * polygon[i][1]
        for i in range(len(polygon) - 1)
    )

def polygon_contains_point(polygon, point):
    # Adapted from
    # https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
    x, y = point
    
    contained = False
    for i in range(1, len(polygon)):
        xi, yi = polygon[i]
        xj, yj = polygon[i-1]
        if ((yi > y) != (yj > y)) and (
            x < xi + (xj - xi) * (y - yi) / (yj - yi)
        ):
            contained = not contained
    return contained

def d2r(p):
    return p[0] * pi / 180.0, p[1] * pi / 180.0

def spherical_distance(p1, p2, r = 6378137.0):
    try:
        londiff = p1[0] - p2[0]
        latdiff = p1[1] - p2[1]

        return r * 2.0 * asin(sqrt(
            sin(latdiff / 2.0) ** 2.0 + cos(p1[1]) * cos(p2[1]) * sin(londiff / 2.0) ** 2.0
        ))
    except Exception as e:
        print(londiff, latdiff, p1, p2)
        raise e

def euclidean_squared_distance(p1, p2):
    xdiff = p1[0] - p2[0]
    ydiff = p1[1] - p2[1]
    
    return xdiff * xdiff + ydiff * ydiff

def euclidean_distance(p1, p2):
    xdiff = p1[0] - p2[0]
    ydiff = p1[1] - p2[1]
    
    return sqrt(xdiff * xdiff + ydiff * ydiff)

# Silly Web Mercator nonsense
# Adapted from
# https://alastaira.wordpress.com/2011/01/23/the-google-maps-bing-maps-spherical-mercator-projection/

def wgs84_to_espg3857(lon, lat):
    x = lon * 20037508.34 / 180;
    y = 20037508.34 / 180 * log(
      tan((90.0 + lat) * pi / 360.0)
    ) / (pi / 180.0)
    
    return x, y

def espg3857_to_wgs84(x, y):
    lon = (x / 20037508.34) * 180
    lat = (y / 20037508.34) * 180
    lat_adj = 180.0 / pi * (2 * atan(exp(lat * pi / 180.0)) - pi / 2.0)
    
    return lon, lat_adj
