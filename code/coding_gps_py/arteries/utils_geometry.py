import math
import numpy as np
from geopy.distance import vincenty
# import gmplot
from pprint import pprint

# average horiz/vertical
COORD2METERS_SCALE = (108811.8 + 110622.9) / 2.0

# Given a line with coordinates 'start' and 'end' and the
# coordinates of a point 'pnt' the proc returns the shortest
# distance from pnt to the line and the coordinates of the
# nearest point on the line.
#
# 1  Convert the line segment to a vector ('line_vec').
# 2  Create a vector connecting start to pnt ('pnt_vec').
# 3  Find the length of the line vector ('line_len').
# 4  Convert line_vec to a unit vector ('line_unitvec').
# 5  Scale pnt_vec by line_len ('pnt_vec_scaled').
# 6  Get the dot product of line_unitvec and pnt_vec_scaled ('t').
# 7  Ensure t is in the range 0 to 1.
# 8  Use t to get the nearest location on the line to the end
#    of vector pnt_vec_scaled ('nearest').
# 9  Calculate the distance from nearest to pnt_vec_scaled.
# 10 Translate nearest back to the start/end line.
# Malcolm Kesson 16 Dec 2012


def vec_dot(v, w):
    x, y = v
    X, Y = w
    return x * X + y * Y


def vec_length(v):
    x, y = v
    return math.sqrt(x * x + y * y)


def vec_vector(b, e):
    x, y = b
    X, Y = e
    return (X - x, Y - y)


def vec_unit(v):
    x, y = v
    mag = vec_length(v)
    return (x / mag, y / mag)


def vec_distance(p0, p1):
    return vec_length(vec_vector(p0, p1))


def vec_scale(v, sc):
    x, y = v
    return (x * sc, y * sc)


def vec_add(v, w):
    x, y = v
    X, Y = w
    return (x + X, y + Y)


def pnt2line(pnt, start, end):
    line_vec = vec_vector(start, end)
    pnt_vec = vec_vector(start, pnt)
    line_len = vec_length(line_vec)

    # if the line is degenerate
    if line_len * COORD2METERS_SCALE < 20:
        return vec_length(pnt_vec), start

    line_unitvec = vec_unit(line_vec)
    pnt_vec_scaled = vec_scale(pnt_vec, 1.0 / line_len)
    t = vec_dot(line_unitvec, pnt_vec_scaled)
    if t < 0.0:
        t = 0.0
    elif t > 1.0:
        t = 1.0
    nearest = vec_scale(line_vec, t)
    dist = vec_distance(nearest, pnt_vec)
    nearest = vec_add(nearest, start)
    return dist, nearest

# EOF # Malcolm Kesson 16 Dec 2012


def dist_pt2route(point, route):
    # distance in meters from point to route
    # using lat long
    # return in meters

    segments = [(route[idx], route[idx + 1]) for idx in range(len(route) - 1)]
    distances = [pnt2line(point, start, end)[0] * COORD2METERS_SCALE for start, end in segments]
    return min(distances)


def fine_route(route, max_dist=100, snap_to_roads=False):
    """
    Refine a route such that no edge is longer than 100m
     to do this we keep adding points along the edges
    """
    # using lat long

    if snap_to_roads:
        # snap to roads google api
        pass

    route_fine = []
    for idx, point in enumerate(route[1:]):
        point_prev = route[idx]
        dist2prev = vincenty(point, point_prev).meters
        if dist2prev > max_dist:
            n_max = max(1, math.ceil(dist2prev / max_dist))
            points_to_add = [tuple((1 - i / n_max) * np.array(point_prev) + i / n_max * np.array(point))
                             for i in range(n_max)]
            route_fine += points_to_add
        else:
            route_fine.append(point_prev)
    route_fine.append(route[-1])
    return route_fine


def trips_from_points(points):
    trips = []
    trip = []
    for idx, point in enumerate(points[:-1]):
        after_prev = points[idx+1]
        dist2after = vincenty(point, after_prev).meters
        if dist2after <= 100:
            trip.append(point)
        else:
            trip.append(point)
            trips.append(trip)
            trip = []

    trips = [trip for trip in trips if len(trip)>1]

    return trips


def frac_overlap(trip_route, network, eps_dist=None):
    """
    This function computes the fraction of ***trip_route*** that is within eps_dist of ***network***
         In other words, for each point on route1, we find if the distance to route2 is > or < eps_dist
         This function returns the fraction of points with distance to route2 < eps_dist

    :param trip_route: list of (lat,long) pairs - these are the points on route1
    :param route2: list of (lat,long) pairs - these are the points on route1
    :param eps_dist: tolerance in meters. Otherwise if None the tolerance is automatically computed.
    :return:
    """
    # using lat long
    # eps_dist in meters

    # overlap tolerance from =100 to =1000 from 0km to 20km
    # only if the parameter eps_dist is not passed as argument (otherwise this is skipped)
    if eps_dist is None:
        total_route_dist = sum([vincenty(trip_route[i], trip_route[i + 1]).meters for i in range(len(trip_route) - 1)])
        eps_dist = max(100, min(1000, total_route_dist / 20))
        print("OVERLAP - using eps = " + str(eps_dist) + " for route length = " + str(total_route_dist))

    # add points along trip_route along the network so that we have points every 100 meters at least
    trip_route_fine = fine_route(trip_route)

    # for each point in (fine) route 1, compute the distances to every edge and get the minimum distance
    overlap_trips = []
    non_overlap_trips = []
    min_distances = []
    overlap_ongoing = True
    curr_trip = []
    for point in trip_route_fine:
        distances = [dist_pt2route(point, edge['edge_route']) for edge in network]
        min_dist = min(distances)
        min_distances.append(min_dist)
        if overlap_ongoing:
            if min_dist <= eps_dist:
                curr_trip.append(point)
            else:
                # curr_trip.append(point)
                overlap_trips.append(curr_trip)
                # todo: add previous point?
                curr_trip = [point]
                overlap_ongoing = False
        else:
            if min_dist > eps_dist:
                curr_trip.append(point)
            else:
                # curr_trip.append(point)
                non_overlap_trips.append(curr_trip)
                # todo: add previous point?
                curr_trip = [point]
                overlap_ongoing = True

    # add last point
    if overlap_ongoing:
        overlap_trips.append(curr_trip)
    else:
        non_overlap_trips.append(curr_trip)

    # for each point in (fine) route 1, 1 if it's close to any edge, 0 if not (threshold = eps_dist)
    close_to_edges = [int(md < eps_dist) for md in min_distances]

    # todo: total distances of overlap and non-overlap
    distances = [vincenty(trip_route_fine[idx], trip_route_fine[idx+1]) for idx in range(len(trip_route_fine) - 1)]
    distances.append(vincenty(trip_route_fine[-1], trip_route_fine[-2]))
    
    # compute fraction of points on trip_route that are close to any edge
    frac_olap_1 = np.mean(close_to_edges)

    return frac_olap_1, overlap_trips, non_overlap_trips


def fine_route_exactly(route, max_dist=300, max_n=50):
    route_fine = fine_route(route)
    point_prev = route_fine[0]
    new_route_fine = [point_prev]
    distance = 0
    for point in route_fine:
        distance += vincenty(point_prev, point).meters
        if distance > max_dist:
            new_route_fine.append(point)
            # reset distance
            distance = 0
        point_prev = point

    if len(new_route_fine) > max_n:
        # print("too many points: " + str(len(new_route_fine)))
        route_fine = fine_route(route, max_dist=50)
        every_n = max(int(math.floor(len(route_fine) / max_n)), 1)
        # print("picking every " + str(every_n))
        # sieve
        new_route_fine = route_fine[0::every_n]

    return new_route_fine



