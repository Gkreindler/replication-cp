import collections
import csv
import datetime
import math
import numpy as np
import statistics
from geopy.distance import great_circle
np.seterr(all='raise')

# assume all times are Indian Standard Time
TIMEZONE = datetime.timezone(datetime.timedelta(hours=5, minutes=30), 'IST')

# in meters
EARTH_RADIUS = 6371010

# text encoding of csv files
ENCODING = 'iso-8859-1'


class Point(object):
    def __init__(self, deviceid, longitude, latitude, accuracy_level=-1, event_time=None, event_year=2016,
                 event_month=1, event_day=1, event_hour=0, event_min=0, event_sec=0, location_code=-1,
                 battery_status=-1, uidp='', name='', redundant=False, id="", idx=-1, **_):
        self.deviceid = str(deviceid)
        self.uidp = str(uidp)
        self.name = str(name)
        self.longitude = float(longitude)
        self.latitude = float(latitude)
        self.accuracy_level = int(accuracy_level)
        self.battery_status = int(battery_status)
        # self.timestamp = datetime.datetime(int(event_year), int(event_month),
        #                                    int(event_day), int(event_hour), int(event_min),
        #                                    int(event_sec), tzinfo=TIMEZONE)
        # todo: switch to read from event_time
        self.timestamp = datetime.datetime.strptime(event_time[2:], '%Y-%m-%d %H:%M:%S')

        self.location_code = int(location_code)  # pre-coded location ID (eg home, work). -1 = none
        self.trip = False  # point is part of a trip (0,1)
        self.redundant = redundant
        self.low_accuracy = (self.accuracy_level > 200)
        self.pid = str(id)
        self.idx = int(idx)

        self.drop = False
        self.drop_only_color = False
        self.drop_type = ''

        self.score_after = -1
        self.score_before = -1

        self.ignore = False

        self.time2prev = -1  # seconds
        self.dist2prev = -1.0
        self.time2next = -1  # seconds
        self.dist2next = -1.0
        self.angle = -1.0  # angles?

        self.markercolor = ""
        self.linecolor = ""

    def __hash__(self):
        return hash((self.longitude, self.latitude, self.accuracy_level,
                     self.timestamp))

    def __repr__(self):
        return '{}({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})'.format(
            self.__class__.__name__, self.deviceid, self.longitude, self.latitude,
            self.accuracy_level, *self.timestamp.timetuple(), self.location_code)

    def __str__(self):
        return '{0.timestamp}: {0.latitude}, {0.longitude}'.format(self)

    def __le__(self, other):
        return self.timestamp <= other.timestamp

    def __lt__(self, other):
        return self.timestamp < other.timestamp


class Location(Point):
    def __init__(self, deviceid, points, start_time, end_time, longitude, latitude, accuracy_level=-1, event_year=2016,
                 event_month=1, event_day=1, event_hour=0, event_min=0, event_sec=0, location_code=-1, **_):

        super(Location, self).__init__(deviceid, longitude, latitude, accuracy_level, event_year,
                                       event_month, event_day, event_hour, event_min, event_sec, location_code)
        self.points = points  # list of points
        self.start_time = start_time
        self.end_time = end_time

    def updateloc(self, devid_mainlocs_dict, min_lc_frequency=0.6):
        """
        Update location code, lat, long and accuracy
        :param devid_mainlocs_dict: Ensure that mainlocs_dict[devid] is passed
        :param min_lc_frequency:
        :return:
        """

        self.updateloc_code()

        # if unknown location code (-1)
        if self.location_code == -1:
            # update the lat and long of the location, using the median of accurate points
            acc_points = [point for point in self.points if point.accuracy_level < 200 and not point.drop]
            if acc_points:
                self.latitude = statistics.median([point.latitude for point in acc_points])
                self.longitude = statistics.median([point.longitude for point in acc_points])

                # all locations have accuracy = 200
                self.accuracy_level = 200

                # median accuracy or at least 200
                # self.accuracy_level = max(200, statistics.median([point.accuracy_level for point in acc_points]))

                # accuracy should be the max radius from the next center, but no less than 100 meters
                # no location should have an accuracy > 1km
                # self.accuracy_level = max([distance(self, point) for point in acc_points])
                # self.accuracy_level = min(max(self.accuracy_level, 150), 1000)

        # known location code
        else:
            # use the values for that location code
            self.latitude = devid_mainlocs_dict[self.location_code].latitude
            self.longitude = devid_mainlocs_dict[self.location_code].longitude
            self.accuracy_level = devid_mainlocs_dict[self.location_code].accuracy_level

    def merge_locations(self, other):
        """
        Take two consecutive locations and create a single location  --- WORK IN PROGRESS
        """
        # Check same location
        assert not distinct(self, other)
        merged_location = self
        merged_location.points += other.points
        return merged_location

    def updateloc_code(self):
        # update the location_code using the following rule: if there is a plurarity of >=60% of a known location code,
        #  then pick that. Otherwise, -1 (default)

        loc_code_freq_dict = {}
        for point in self.points:
            if point.location_code in loc_code_freq_dict:
                loc_code_freq_dict[point.location_code] += 1
            else:
                loc_code_freq_dict[point.location_code] = 1

        possible_code = -1
        max_freq = 0

        for code, freq in loc_code_freq_dict.items():
            if freq > max_freq:
                possible_code = code
                max_freq = freq

        if max_freq/float(len(self.points)) >= 0.6:
            self.location_code = possible_code
        else:
            self.location_code = -1


class Trip(object):
    """
     Contains information on a trip.
     points is a list of points strictly in the trip (not including the origin and destination)
     orig is the starting location with the timestamp of the last known location at that location (before trip)
     dest is the ending location with the timestamp of the first known location at that location (after trip)
    """
    def __init__(self, deviceid, orig, dest, points, start_time=-1, start_precision=-1,
                 end_time=-1, end_precision=-1,name='', uidp='', charge=-1, **_):
        self.name = name
        self.uidp = uidp
        self.deviceid = deviceid
        self.orig = orig # location before this trip started
        self.dest = dest # location before this trip started
        self.points = points  # list of points in trip.
        self.charge = charge # charge for this trip in the treatment

        self.start_time = start_time
        self.start_precision = start_precision
        self.end_time = end_time
        self.end_precision = end_precision

    def starttime(self):
        """
        Compute estimated start time as midpoint between origin and first point timestamps
        :return: estimated start time, and precision (difference between origin and first point)
        """
        try:
            startdelta = self.points[0].timestamp - self.orig.timestamp
            lambda_orig = 0.5
            return [self.orig.timestamp + lambda_orig * startdelta, startdelta.total_seconds()]
        except:
            return [self.start_time, self.start_precision]

    def endtime(self):
        """
        Compute estimated end time as midpoint between last point and destination timestamps
        :return: estimated end time, and precision (difference between last point and destination)
        """
        try:
            enddelta = self.dest.timestamp - self.points[-1].timestamp
            lambda_dest = 0.5
            return [self.points[-1].timestamp + (1-lambda_dest) * enddelta, enddelta.total_seconds()]
        except:
            return [self.end_time, self.end_precision]

    def date(self):
        return self.starttime()[0].date()

    def od(self, oid, did):
        return self.orig.location_code == oid and self.dest.location_code == did

    def total_meters(self):
        total_dist = 0
        try:
            prev_point = self.orig
            for point in self.points:
                total_dist += min(500, distance(prev_point, point))
                prev_point = point
            total_dist += distance(prev_point, self.dest)
            return total_dist
        except Exception as e:
            print("Error while calculating total trip distancei in meters:")
            print(e)
            return -1


def distance(a, b):
    return great_circle((a.latitude, a.longitude), (b.latitude, b.longitude)).meters


def distance_dict(a, b):
    return great_circle((float(a['latitude']), float(a['longitude'])),
                        (float(b['latitude']), float(b['longitude']))).meters


def distinct(self, other, tol=1.1):
    if distance(self, other) > tol * max(self.accuracy_level, other.accuracy_level):
        return True
    else:
        return False


def dist_point_array(point_source, point_array):
    """Calculate the distance in meters from one Point to several points.

    This uses the equirectangular approximation, so it only works if the points
    are relatively close together.
    """

    latitudes = np.array([pt.latitude for pt in point_array])
    longitudes = np.array([pt.longitude for pt in point_array])

    dy = np.radians(point_source.latitude - latitudes)
    cosy = math.cos(point_source.latitude)
    dx = np.radians(point_source.longitude - longitudes) * cosy

    return np.sqrt(dx ** 2 + dy ** 2) * EARTH_RADIUS


def angle(p1, p2, p3):
    """
    Angle at p2 in degrees
    :param p1:
    :param p2:
    :param p3:
    :return:
    """
    p12 = distance(p1, p2)
    p13 = distance(p1, p3)
    p23 = distance(p2, p3)

    if p12 * p23 == 0:
        return -999.0
    try:
        temp = (p12 ** 2 + p23 ** 2 - p13 ** 2) / (2 * p12 * p23)
        # assert temp < 1.000000001
        # assert temp > -1.000000001
        temp = max(-1, min(1, temp))
        return np.degrees(np.arccos(temp))
    except Exception as e:
        print(p1)
        print(p2)
        print(p3)
        print(temp)
        raise e


def point_stats(points_list, debug=False):
    """
    Compute dist2prev, dist2next, time2prev, time2next, and angle
    IGNORE dropped points!
    :param points_list:
    :return:
    """

    # index of prev/next not missing point
    idx_prev = [-1] * len(points_list)
    idx_next = [-1] * len(points_list)

    idx = 0  # current point
    idx_p1 = -1  # previous
    while idx < len(points_list):
        if not points_list[idx].drop:
            if debug and idx % 100 == 0:
                print('.', end='')
            current_points = [idx]
            idx_n1 = -1

            # keep adding points to "current point" until we reach a different point
            idx2 = idx + 1
            while idx2 < len(points_list):
                if not points_list[idx2].drop:
                    dist = distance(points_list[idx], points_list[idx2])
                    if dist <= 1:
                        # idx2 is essentially at the same location as idx, so add to list
                        current_points.append(idx2)
                    else:
                        # we've reached a different point, so we can stop
                        idx_n1 = idx2
                        break
                idx2 += 1

            # update all current points:
            for idx2 in current_points:
                idx_prev[idx2] = idx_p1
                idx_next[idx2] = idx_n1

            # if a previous point exists
            # if idx_p1 >= 0:
            #     idx_next[idx_p1] = idx

            # last point in current (non-dropped) point list becomes the "previous point"
            idx_p1 = current_points[-1]

            # skip all point in the current location (already covered)
            idx = current_points[-1]

        idx += 1

    idx = 0
    for idx, point in enumerate(points_list):
        if not point.drop:

            # previous point exists
            if idx_prev[idx] >= 0:
                point_prev = points_list[idx_prev[idx]]
                points_list[idx].dist2prev = distance(point, point_prev)
                points_list[idx].time2prev = (point.timestamp - point_prev.timestamp).total_seconds()

            # next point exists
            if idx_next[idx] >= 0:
                point_next = points_list[idx_next[idx]]
                points_list[idx].dist2next = distance(point, point_next)
                points_list[idx].time2next = (point_next.timestamp - point.timestamp).total_seconds()

            # both previous and next exist
            if idx_prev[idx] >= 0 and idx_next[idx] >= 0:
                point_prev = points_list[idx_prev[idx]]
                point_next = points_list[idx_next[idx]]
                points_list[idx].angle = angle(point_prev, point, point_next)
    # print('!')


def dict2point(my_dict):

    hh = int(my_dict['time'][0:2])
    mm = int(my_dict['time'][3:5])
    ss = int(my_dict['time'][6:8])

    year = int(my_dict['date'][0:4])
    month = int(my_dict['date'][5:7])
    day = int(my_dict['date'][8:10])

    pt = Point(deviceid=my_dict.get("deviceid", ""),
               latitude=float(my_dict.get("latitude", 0)),
               longitude=float(my_dict.get("longitude", 0)),
               accuracy_level=int(my_dict.get("accuracy_level", 200)),
               event_year=year,
               event_month=month,
               event_day=day,
               event_hour=hh,
               event_min=mm,
               event_sec=ss)
    return pt


def generate_points(filename):
    """ Generate a list of Points from a csv file. """
    try:
        with open(filename, encoding=ENCODING) as csvfile:
            yield from (Point(**d) for d in csv.DictReader(csvfile) if abs(float(d.get('latitude', 0.0))) > 0.01)
    except Exception as e:
        raise e
        # return []


def sort_by_date(points, date_min=None, date_max=None, dates=None):
    """Sorts a list of Points by deviceid and date."""

    # no dates -> full range
    if not dates:
        if not date_min:
            date_min = datetime.date(2000, 1, 1)
        if not date_max:
            date_max = datetime.date(2100, 1, 1)
        dates = [date_min + datetime.timedelta(days=i) for i in range((date_max - date_min).days + 1)]

    # separate into dates
    out = {}
    for p in points:
        date = p.timestamp.date()
        if date in dates:
            if date not in out:
                out[date] = []
            out[date].append(p)

    # sort by time within each date
    for date in sorted(out.keys()):
        # uses the pointclass point comparison: point1 < point2 iff timestamp1 < timestamp2
        assert all(out[date][i] <= out[date][i+1] for i in range(len(out[date])-1))
        # assert all(out[date][i].pid <= out[date][i + 1].pid for i in range(len(out[date]) - 1))  # NOT ALWAYS TRUE
        out[date] = sorted(out[date])

        # add idx from 1 to N_POINTS
        for idx, point in enumerate(out[date]):
            point.idx = idx  # todo: +1 or not +1

    return out


def point_list_center_2(point_list, accuracy_fraction=80, min_acc=200):
    from statistics import median
    from geopy.distance import vincenty
    from geopy.distance import great_circle
    from numpy import percentile

    center_lat = median([float(point.latitude) for point in point_list])
    center_lon = median([float(point.longitude) for point in point_list])
    try:
        distances = [vincenty((center_lat, center_lon), (point.latitude, point.longitude)).meters for point in point_list]
    except:
        distances = [great_circle((center_lat, center_lon), (point.latitude, point.longitude)).meters for point in
                     point_list]

    # accuracy that covers % of points
    center_acc = max(min_acc, percentile(distances, accuracy_fraction))

    return center_lat, center_lon, center_acc


def point_list_center(point_list, accuracy_fraction=80, min_acc=200):
    from statistics import median
    # from geopy.distance import vincenty
    from geopy.distance import great_circle
    # from numpy import percentile
    import numpy as np

    center_lat = median([float(point['latitude']) for point in point_list])
    center_lon = median([float(point['longitude']) for point in point_list])
    try:
        distances = [great_circle((center_lat, center_lon), (point['latitude'], point['longitude'])).meters for point in point_list]
    except:
        distances = [great_circle((center_lat, center_lon), (point['latitude'], point['longitude'])).meters for point in
                     point_list]

    # accuracy that covers % of points
    center_acc = max(min_acc, np.percentile(distances, accuracy_fraction))

    # std of all points
    center_std = np.array(distances).std()

    return center_lat, center_lon, center_acc, center_std
