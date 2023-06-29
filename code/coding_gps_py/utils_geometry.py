# convex hull (Graham scan by x-coordinate) and diameter of a set of points
# David Eppstein, UC Irvine, 7 Mar 2002
from geopy.distance import great_circle
import math
from math import cos, pi, pow
TILE_SIZE = 256.0
ZOOM_LEVEL = 15
pixelOrigin = (TILE_SIZE / 2, TILE_SIZE / 2)
pixelPerLonDeg = TILE_SIZE / 360.0
pixelPerLonRad = TILE_SIZE / (2 * math.pi)


def orientation(p, q, r):
    '''Return positive if p-q-r are clockwise, neg if ccw, zero if colinear.'''
    return (q[1] - p[1]) * (r[0] - p[0]) - (q[0] - p[0]) * (r[1] - p[1])


def hulls(Points):
    '''Graham scan to find upper and lower convex hulls of a set of 2d points.'''
    U = []
    L = []
    Points.sort()
    for p in Points:
        while len(U) > 1 and orientation(U[-2], U[-1], p) <= 0: U.pop()
        while len(L) > 1 and orientation(L[-2], L[-1], p) >= 0: L.pop()
        U.append(p)
        L.append(p)
    return U, L


def rotatingCalipers(Points):
    '''Given a list of 2d points, finds all ways of sandwiching the points
between two parallel lines that touch one point each, and yields the sequence
of pairs of points touched by each pair of lines.'''
    U, L = hulls(Points)
    i = 0
    j = len(L) - 1
    while i < len(U) - 1 or j > 0:
        yield U[i], L[j]

        # if all the way through one side of hull, advance the other side
        if i == len(U) - 1:
            j -= 1
        elif j == 0:
            i += 1

        # still points left on both lists, compare slopes of next hull edges
        # being careful to avoid divide-by-zero in slope calculation
        elif (U[i + 1][1] - U[i][1]) * (L[j][0] - L[j - 1][0]) > \
                        (L[j][1] - L[j - 1][1]) * (U[i + 1][0] - U[i][0]):
            i += 1
        else:
            j -= 1


def diameter(Points):
    '''Given a list of 2d points, returns the pair that's farthest apart.'''
    diam, pair = max([((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2, (p, q))
                      for p, q in rotatingCalipers(Points)])

    return great_circle(pair[0], pair[1]).meters, pair


""" formerly from utils_geometry in bang_treat repo """


def bound(val, valMin, valMax):
    res = max(val, valMin)
    res = min(res, valMax)
    return res


def latlon_to_coord(latlon, zoom):
    lat, lon = latlon
    x = pixelOrigin[0] + lon * pixelPerLonDeg

    siny = bound(math.sin(lat * math.pi / 180.0), -0.99999, 0.9999)
    y = pixelOrigin[1] + 0.5 * math.log((1 + siny) / (1 - siny)) * - pixelPerLonRad
    c = math.pow(2, zoom)
    return x * c, y * c


def coord_to_latlon(coord, zoom):
    c = math.pow(2, zoom)
    x, y = coord
    x1 = x / c
    y1 = y / c
    lon = (x1 - pixelOrigin[0]) / pixelPerLonDeg
    latRadians = (y1 - pixelOrigin[1]) / -pixelPerLonRad
    lat = (2 * math.atan(math.exp(latRadians)) - math.pi / 2) / (math.pi / 180)
    return lat, lon


def get_intersection(vec1, vec2):
    # take in vec1=((sx1,sy1),(ex1,ey1)) and vec2=((sx2,xy2),(ex2,ey2))
    # return one intersection (if exist) (ix1,iy1) of
    # the lines form by the above vectors, or None if the two vectors don't intersect

    # special use: if (sx1,sy1)==(ex1,ey1) the function would check whether (sx1,sy1) is on the line form by vec2
    # and if yes return (sx1,sy1)
    s1, e1 = vec1
    s2, e2 = vec2
    sx1, sy1 = s1
    ex1, ey1 = e1
    sx2, sy2 = s2
    ex2, ey2 = e2

    diffy1 = ey1 - sy1
    diffx1 = ex1 - sx1
    diffy2 = ey2 - sy2
    diffx2 = ex2 - sx2

    a = diffx2 * diffy1 - diffy2 * diffx1
    b = (sx1 - sx2) * diffy1 - (sy1 - sy2) * diffx1
    if a == 0:
        if b == 0:
            # two vectors are in the same line, can return any point on that line
            return s1
        else:
            # system of linear equations have no solutions
            return None
    else:
        t2 = b * 1.0 / a
        return (sx2 + t2 * diffx2, sy2 + t2 * diffy2)


def get_projection(point, vec):
    # return the projection of point on the line form by vec
    s, e = vec
    sx, sy = s
    ex, ey = e

    diffx = ex - sx
    diffy = ey - sy
    if diffx == 0 and diffy == 0:
        return s
    else:
        perp_vec = (point, (point[0] + diffy, point[1] - diffx))
        return get_intersection(perp_vec, vec)


def squaredist(start, end):
    # return square of the distance between points start and end
    return (start[0] - end[0]) ** 2 + (start[1] - end[1]) ** 2


def check_inside_latlong(vec, point, r):
    lat = point[0]

    vec = latlon_to_coord(vec[0], ZOOM_LEVEL), latlon_to_coord(vec[1], ZOOM_LEVEL)
    point = latlon_to_coord(point, ZOOM_LEVEL)
    ground_resolution = (cos(lat * pi / 180) * 2 * pi * 6378137) / (256 * pow(2, ZOOM_LEVEL))
    r /= ground_resolution

    return check_inside(vec, point, r)


def check_inside(vec, point, r):
    # check if the line segment vec intersect the circle
    # of radius r centered at point, return True if yes
    proj = get_projection(point, vec)
    square_r = r * r
    if squaredist(point, proj) > square_r:
        return False
    else:
        s, e = vec
        dot0 = ((s[0] < proj[0]) and (e[0] > proj[0])) or ((s[0] > proj[0]) and (e[0] < proj[0]))
        dot1 = ((s[1] < proj[1]) and (e[1] > proj[1])) or ((s[1] > proj[1]) and (e[1] < proj[1]))
        if dot0 or dot1:
            return True
        else:
            if squaredist(point, s) <= square_r or squaredist(point, e) <= square_r:
                return True
            else:
                return False


def check_inside_get_intersection(vec, point, r):
    # check if the line segment vec intersect the circle
    # of radius r centered at point,
    # return None if the segment doesn't intersect the circle
    # else return list of the points where the line segment intersect the circle's circumference
    # if the returned list consist of 2 points, then 1st point in list is closer to starting point of the segment,
    # and 2nd point is closer to the end point of the line segment
    proj = get_projection(point, vec)
    square_r = r * r
    square_proj = squaredist(point, proj)
    if square_proj > square_r:
        return None
    else:
        s, e = vec
        dot0 = ((s[0] < proj[0]) and (e[0] > proj[0])) or ((s[0] > proj[0]) and (e[0] < proj[0]))
        dot1 = ((s[1] < proj[1]) and (e[1] > proj[1])) or ((s[1] > proj[1]) and (e[1] < proj[1]))
        distS = squaredist(point, s)
        distE = squaredist(point, e)
        flag = False
        if dot0 or dot1:
            flag = True
        else:
            if distS <= square_r or distE <= square_r:
                flag = True
        if flag:
            ans = []
            # !!:zero division error]
            square_d = squaredist(s, e)
            if square_d == 0:
                # s and e is the same point inside the circle
                if (squaredist(point, s) == square_r):
                    return [s]
                else:
                    return []
            ratio = math.sqrt((square_r - square_proj) * 1.0 / square_d)
            diffx = s[0] - e[0]
            diffy = s[1] - e[1]
            if distS >= square_r:
                interS = (proj[0] + ratio * diffx, proj[1] + ratio * diffy)
                ans.append(interS)
            if distE >= square_r:
                interE = (proj[0] + ratio * -diffx, proj[1] + ratio * -diffy)
                ans.append(interE)
            return ans
        else:
            return None


""" formerly from utils_geometry_2 in repo bang_treat """

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
        dist = vec_length(pnt_vec) * COORD2METERS_SCALE
        return dist, start

    line_unitvec = vec_unit(line_vec)
    pnt_vec_scaled = vec_scale(pnt_vec, 1.0 / line_len)
    t = vec_dot(line_unitvec, pnt_vec_scaled)
    if t < 0.0:
        t = 0.0
    elif t > 1.0:
        t = 1.0
    nearest = vec_scale(line_vec, t)
    dist = vec_distance(nearest, pnt_vec) * COORD2METERS_SCALE
    nearest = vec_add(nearest, start)
    return dist, nearest

# EOF # Malcolm Kesson 16 Dec 2012
