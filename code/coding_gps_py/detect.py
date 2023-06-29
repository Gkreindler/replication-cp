# detect trips from a list of points
# the list contains points from the same day, ordered by event time

import pointclass as pc
import datetime
from utils import time_color, clamp, str_color
from geopy.distance import great_circle
import collections


def sharp_angle_real_trip(idx_to_drop, points_list):
    if len(idx_to_drop) < 3:
        return False
    p1 = points_list[idx_to_drop[-1]]
    p2 = points_list[idx_to_drop[-2]]
    p3 = points_list[idx_to_drop[-3]]
    if len(idx_to_drop) >= 4:
        p4 = points_list[idx_to_drop[-4]]
    else:
        p4 = p3

    d12 = pc.distance(p1, p2)
    d23 = pc.distance(p2, p3)
    d34 = pc.distance(p3, p4)
    d14 = pc.distance(p1, p4)

    score = 0
    score += 50 <= d12 <= 500
    score += 50 <= d23 <= 500
    score += 50 <= d34 <= 500
    score += 2 * (d12 + d23 + d34 <= 1.5 * d14)

    return score >= 3


def above_speed_limit(dist, dur, conservative=True):
    """ distance in meters, duration in seconds """
    assert dist >= 0
    assert dur > 0
    if conservative:
        if dist < 1000:
            # 50m/s = 180 km/h
            speed_threshold = 50
            return dist / dur > speed_threshold
        elif dist < 2000:
            # 50m/s = 180 km/h
            # 33m/s = 120 km/h
            speed_threshold = ((dist - 1000) * 33 + (2000 - dist) * 50) / 1000
            assert 33 <= speed_threshold <= 50
            return dist / dur > speed_threshold
        elif dist < 4000:
            # 33m/s = 120 km/h
            # 22m/s = 80 km/h
            speed_threshold = ((dist - 2000) * 22 + (4000 - dist) * 33) / 2000
            assert 22 <= speed_threshold <= 33
            return dist / dur > speed_threshold
        else:
            # 22m/s = 80 km/h
            speed_threshold = 22
            return dist / dur > speed_threshold
    else:
        if dist < 1000:
            # 8.3m/s = 30 km/h
            speed_threshold = 8.3
            return dist / dur > speed_threshold
        elif dist < 2000:
            # 8.3m/s = 30 km/h
            # 6.94m/s = 25 km/h
            speed_threshold = ((dist - 1000) * 6.94 + (2000 - dist) * 8.3) / 1000
            assert 6.94 <= speed_threshold <= 8.3
            return dist / dur > speed_threshold
        elif dist < 4000:
            # 6.94m/s = 25 km/h
            # 5.55m/s = 20 km/h
            speed_threshold = ((dist - 2000) * 5.55 + (4000 - dist) * 6.94) / 2000
            assert 5.55 <= speed_threshold <= 6.94
            return dist / dur > speed_threshold
        elif dist < 8000:
            # 5.55m/s = 20 km/h
            # 4.16m/s = 15 km/h
            speed_threshold = ((dist - 4000) * 4.16 + (8000 - dist) * 5.55) / 4000
            assert 4.16 <= speed_threshold <= 5.55
            return dist / dur > speed_threshold
        else:
            # 4.16m/s = 15 km/h
            speed_threshold = 4.16
            return dist / dur > speed_threshold


def get_prev_idx(idx, points_list):
    if not 0 <= idx < len(points_list):
        return -1

    idx_new = idx
    while idx_new > 0:
        idx_new -= 1
        if not points_list[idx_new].drop:
            return idx_new

    # it means all previous points were drops
    return idx


def get_next_idx(idx, points_list):
    if not 0 <= idx < len(points_list):
        return -1

    idx_new = idx
    while idx_new < len(points_list) - 1:
        idx_new += 1
        if not points_list[idx_new].drop:
            return idx_new
    # it means all subsequent points were drops
    return idx


def main_locations(allpoints, allmainlocs):
    """
    Classify points based on whether they are one of the main pre-coded main locations
    Include location score for each point
    Include distance and time to previous
    :param allpoints:
    :param allmainlocs:
    :return:
    """
    for devid in allpoints.keys():
        if devid in allmainlocs:  # maybe find more computationally effective way of doing this
            points = allpoints[devid]
            mainlocs = allmainlocs[devid]
            for date in sorted(points.keys()):
                points_today = points[date]

                for idx, point in enumerate(points_today):
                    # point.dist2prev = -1
                    # point.time2prev = -1

                    # if idx > 0:
                        # allpoints[devid][date][idx].dist2prev = pc.distance(point, points_today[idx - 1])
                        # timeinsec = (point.timestamp - points_today[idx - 1].timestamp).total_seconds()
                        # allpoints[devid][date][idx].time2prev = math.floor(timeinsec/6)/10
                    for loc_id in mainlocs:
                        mainloc = mainlocs[loc_id]
                        if pc.distance(point, mainloc) < 1.2 * max(mainloc.accuracy_level, point.accuracy_level):
                            allpoints[devid][date][idx].location_code = loc_id
            print('.', end="")
    print('\n')
    return allpoints


def dist_ramp(p1, p2, d1=100, d2=200):
    dist = pc.distance(p1, p2)
    return (dist <= d1) + (d1 < dist < d2) * (d2 - dist) / (d2 - d1)


def kernel_weight(time_dist, dt0=1.5, dt1=5.0, dt2=15.0):
    return (time_dist - dt0) / (dt1 - dt0) * (dt0 <= time_dist < dt1) + \
           (dt2 - time_dist) / (dt2 - dt1) * (dt1 <= time_dist < dt2)


def kernel_score(point, pointlist, dt0=1.5, dt1=5.0, dt2=15.0):

    # time to point (in minutes)
    time_dists = [abs((pt.timestamp - point.timestamp).total_seconds() / 60) for pt in pointlist]
    # weights depend on the time to point
    weights = [kernel_weight(time_dist, dt0=dt0, dt1=dt1, dt2=dt2) for time_dist in time_dists]
    tot_weight = sum(weights)

    # score(distance(pt,point))
    dist_scores = [dist_ramp(point, pt) for pt in pointlist]

    if tot_weight:
        my_kernel_score = sum([i*j for (i, j) in zip(weights, dist_scores)]) / tot_weight
        return my_kernel_score, tot_weight
    else:
        return 0.5, tot_weight


def find_points_after(idx, points_list, after=True, dt0=1.5, dt1=5.0, dt2=15.0):
    # if after:
    idx_list = range(idx + 1, len(points_list))
    # else:
    #     idx_list = range(0, idx)

    point = points_list[idx]

    pts_before_0 = [i for i in idx_list if
                    abs((points_list[i].timestamp - point.timestamp).total_seconds()) < dt0 * 60]

    pts_before_1 = [i for i in idx_list if
                    dt0 * 60 <= abs((points_list[i].timestamp - point.timestamp).total_seconds()) < dt1 * 60]

    pts_before_2 = [i for i in idx_list if
                    dt1 * 60 <= abs((points_list[i].timestamp - point.timestamp).total_seconds()) < dt2 * 60]

    pts_before_3 = [i for i in idx_list if
                    dt2 * 60 <= abs((points_list[i].timestamp - point.timestamp).total_seconds())]

    if pts_before_2 or not pts_before_3:

        # apply kernel to 1.5 - 15 if non-empty
        if pts_before_1 + pts_before_2:
            my_list = pts_before_1 + pts_before_2
            return my_list[-1]
            # return my_list[0]

        # or try only before 1.5 (if we're here we know pts_before_3 is empty, so this is all)
        elif pts_before_0:
            return pts_before_0[-1]
            # return pts_before_0[-1]

        # maybe you're the first point
        else:
            return []

    # 1.5 - 15 empty but there's something beyond 15
    else:
        # score with first far-away points
        # if after:
        return pts_before_3[0]
        # else:
        # return pts_before_3[-1]


def loc_score(idx, points_list, after=True, dt0=1.5, dt1=5.0, dt2=15.0, dt3=30.0):
    """
    segment 0 = 0-1.5
    segment 1 = 1.5-5
    segment 2 = 5-15
    segment 3 = 15-30

    :param idx:
    :param points_list:
    :param after:
    :param dt0:
    :param dt1:
    :param dt2:
    :return:
    """
    fpt = None
    point = points_list[idx]

    """ points_list is already sorted by time, so just go up or down """
    pts_0 = []
    pts_1 = []
    pts_2 = []
    pts_3 = []
    if after:
        idx_candidate = get_next_idx(idx, points_list)
    else:
        idx_candidate = get_prev_idx(idx, points_list)

    stagnant = idx == idx_candidate
    if not stagnant:
        fpt = points_list[idx_candidate]

    while abs((points_list[idx_candidate].timestamp - point.timestamp).total_seconds()) < dt0 * 60 \
            and not stagnant:
        pts_0.append(points_list[idx_candidate])

        if after: idx_next = get_next_idx(idx_candidate, points_list)
        else:     idx_next = get_prev_idx(idx_candidate, points_list)

        if idx_next == idx_candidate: stagnant = True
        else: idx_candidate = idx_next

    while dt0 * 60 <= abs((points_list[idx_candidate].timestamp - point.timestamp).total_seconds()) < dt1 * 60 \
            and not stagnant:
        pts_1.append(points_list[idx_candidate])

        if after: idx_next = get_next_idx(idx_candidate, points_list)
        else:     idx_next = get_prev_idx(idx_candidate, points_list)

        if idx_next == idx_candidate: stagnant = True
        else: idx_candidate = idx_next

    while dt1 * 60 <= abs((points_list[idx_candidate].timestamp - point.timestamp).total_seconds()) < dt2 * 60 \
            and not stagnant:
        pts_2.append(points_list[idx_candidate])

        if after: idx_next = get_next_idx(idx_candidate, points_list)
        else:     idx_next = get_prev_idx(idx_candidate, points_list)

        if idx_next == idx_candidate: stagnant = True
        else: idx_candidate = idx_next

    while dt2 * 60 <= abs((points_list[idx_candidate].timestamp - point.timestamp).total_seconds()) < dt3 * 60 \
            and not stagnant:
        pts_3.append(points_list[idx_candidate])

        if after: idx_next = get_next_idx(idx_candidate, points_list)
        else:     idx_next = get_prev_idx(idx_candidate, points_list)

        if idx_next == idx_candidate: stagnant = True
        else: idx_candidate = idx_next

    # fpt = None
    # if after:
    #     pts_list_full = points_list[(idx + 1):]
    #     if pts_list_full:
    #         fpt = pts_list_full[0]
    # else:
    #     pts_list_full = points_list[:idx]
    #     if pts_list_full:
    #         fpt = pts_list_full[-1]
    #
    # point = points_list[idx]
    #
    # pts_0 = [pt for pt in pts_list_full if
    #               abs((pt.timestamp - point.timestamp).total_seconds()) < dt0 * 60]
    #
    # pts_1 = [pt for pt in pts_list_full if
    #               dt0 * 60 <= abs((pt.timestamp - point.timestamp).total_seconds()) < dt1 * 60]
    #
    # pts_2 = [pt for pt in pts_list_full if
    #               dt1 * 60 <= abs((pt.timestamp - point.timestamp).total_seconds()) < dt2 * 60]
    #
    # pts_3 = [pt for pt in pts_list_full if
    #               dt2 * 60 <= abs((pt.timestamp - point.timestamp).total_seconds()) < dt3 * 60]

    # if we have 5-15
    if pts_2:
        # apply kernel to 1.5 - 15 if non-empty
        # if pts_before_1 + pts_before_2:
        k_score, k_tot_w = kernel_score(point, pts_1 + pts_2, dt0=dt0, dt1=dt1, dt2=dt2)
        return k_score

        # or try only before 1.5 (if we're here we know pts_before_3 is empty, so this is ALL)
        # elif pts_before_0:
        #     return kernel_score(point, pts_before_0, dt0=0, dt1=dt0, dt2=dt1)
        
        # maybe you're the first point
        # else:
        #     return 0.5

    # else widen before/after
    else:
        k_score, k_tot_w = kernel_score(point, pts_0 + pts_1 + pts_3, dt0=0, dt1=dt2, dt2=dt3)
        if k_tot_w or not fpt:
            return k_score
        else:
            dist_score = dist_ramp(point, fpt)
            if after:
                return min(dist_score, k_score)
            else:
                return max(dist_score, k_score)

    # 5 - 15 empty but there's something 15-30
    # else:
    #     # score with first far-away points
    #     if after:
    #         fpt = pts_before_3[0]
    #     else:
    #         fpt = pts_before_3[-1]
    #     return kernel_score(point, [fpt], dt0=0, dt1=dt2, dt2=dt3)


def location_score(points_list, dt0=1.5, dt1=5.0, dt2=15.0, d1=100, d2=200):
    for point in points_list:
        point.score_after = -1

    # work with non-dropped points
    points_list_active = [point for point in points_list if not point.drop]

    for i, point in enumerate(points_list_active):

        point.score_after = loc_score(i, points_list_active, after=True, dt0=dt0, dt1=dt1, dt2=dt2)
        point.score_before = loc_score(i, points_list_active, after=False, dt0=dt0, dt1=dt1, dt2=dt2)

    return points_list


def detect_highradius(points_list, min_accuracy_radius=249, debug=False):
    """
    Mark points with high radius for deletion
    Compile some statistics on dropped points
    :param points_list:
    :param min_accuracy_radius:
    :return:
    """
    drop_streak_time_0 = datetime.datetime(year=2017, month=1, day=1, hour=0, minute=0, second=0)
    i_streak = 0
    n_total_drop = 0
    n_streak_drop = 0
    streak_drop = False
    n_drop = 0

    for idx, point in enumerate(points_list):
        if point.accuracy_level > min_accuracy_radius:
            point.drop = True
            point.drop_type = 'low_acc'
            n_drop += 1

            if not streak_drop:
                idx_1 = max(0, idx - 1)
                drop_streak_time_0 = points_list[idx_1].timestamp
                n_streak_drop = 0
                streak_drop = True

            n_streak_drop += 1
            n_total_drop += 1

        else:
            point.drop = False
            if streak_drop:
                i_streak += 1
                drop_streak_time = (point.timestamp - drop_streak_time_0).total_seconds() / 60
                if drop_streak_time > 5:
                    if debug:
                        print("Streak number: " + str(i_streak) + " of " + str(n_streak_drop) + " dropped points.")
                        print("Streak start time: " + str(drop_streak_time_0) + ", total duration (minutes) " + "{0:.2f}".format(drop_streak_time))
                    # print("\n\n")
                streak_drop = False

    print("Dropped " + str(n_drop) + " points out of " + str(len(points_list)))


def detect_crowded_outliers(points_list, min_close_n=2, min_close_dist=5, max_far_n=0, max_far_dist=25, debug=False):
    """ Non-consecutive clusters of points = 2 """

    point_seen = [False] * len(points_list)

    for idx, point in enumerate(points_list):
        # mark all already dropped points as "seen"
        if point.drop:
            point_seen[idx] = True

        # if the angle is high, cannot be a crowded outlier (?)
        if point.angle > 120:
            point_seen[idx] = True

    # go through (unseen) points and check if they are "crowded outliers"
    for idx, point in enumerate(points_list):
        if not point_seen[idx]:
            # point_array = [pt for idx2, pt in enumerate(points_list) if not point_seen[idx2]]
            distances = pc.dist_point_array(point, points_list)
            dist_close = [distance < min_close_dist for distance in distances]

            # count the number of non-consecutive clusters of points that are close to (point)
            # dist_close_n = sum(dist_close)
            dist_close_n = 0
            pt_close = False
            for idx2 in range(len(dist_close)):
                if dist_close[idx2]:
                    # only for first point
                    if not pt_close:
                        dist_close_n += 1
                        pt_close = True
                else:
                    pt_close = False

            dist_far_n = sum([min_close_dist <= distance < max_far_dist for distance in distances])

            if dist_close_n >= min_close_n and dist_far_n <= max_far_n:

                # mark as crowded outlier! (together with all close points)
                for idx2, pt in enumerate(points_list):
                    if dist_close[idx2]:
                        points_list[idx2].drop = True
                        points_list[idx2].drop_type = 'crowded_outlier'
                        points_list[idx2].drop_only_color = False
                        points_list[idx2].markercolor = str_color(0, 255, 0)
                        point_seen[idx2] = True

            point_seen[idx] = True


def detect_sharp_angle(points_list, max_angle=30, min_dist_1=100, min_dist_2=250):
    """
    Detect a sharpt jump. Minimum 100 meters, and the angle threshold (should be below) is
    linear between 100 and 250 meters (between 0 and 30*, and then at 30*.
    :param points_list:
    :param max_angle:
    :param min_dist_1:
    :param min_dist_2:
    :return:
    """

    for idx, point in enumerate(points_list):
        if not point.drop:
            min_dist_prev_next = min(point.dist2prev, point.dist2next)
            if min_dist_prev_next > min_dist_1:

                # linear maximum angle
                if min_dist_prev_next > min_dist_2:
                    max_angle_threshold = max_angle
                else:
                    max_angle_threshold = max_angle * (min_dist_prev_next - min_dist_1) / (min_dist_2 - min_dist_1)
                    assert 0 <= max_angle_threshold <= max_angle

                if point.angle < max_angle_threshold:
                    point.drop = True
                    point.drop_type = 'sharp_angle'
                    point.drop_only_color = False
                    point.markercolor = str_color(255, 0, 0)


def detect_sharp_multiple(points_list, dist_min=750, max_dur=300, go_fwd=True,
                          debug=False, prefix=''):
    """ Detect a sharp angle with multiple points in the destination.
        //speed_min_0 = 33m/s = 120 Km/h
        speed_min_0 = 8.33m/s = 30 Km/h
        speed_min_1 = 22m/s = 90 Km/h
    """

    if not go_fwd:
        points_list = points_list[::-1]

    pts_sharp = []
    for idx, point in enumerate(points_list):
        if not point.drop:
            idx_prev = get_prev_idx(idx, points_list)
            point_prev = points_list[idx_prev]
            dist2prev = pc.distance(point, point_prev)
            time2prev = abs((point.timestamp - point_prev.timestamp).total_seconds())

            # jump condition:
            #  - super short time and distance > 750
            #  - medium time and high speed and and distance > 750
            #  - long time and very high speed
            is_jump = False
            if time2prev <= 10 and dist2prev > dist_min: is_jump = True
            if time2prev > 10 and above_speed_limit(dist2prev, time2prev,
                                                    conservative=False) and dist2prev > dist_min: is_jump = True
            if time2prev > 10 and above_speed_limit(dist2prev, time2prev,
                                                    conservative=True) and dist2prev > 0.5 * dist_min: is_jump = True

            if is_jump:
                # find if there is a return to point_prev in next max_dur seconds
                comes_back = False
                dur = 0
                idx2 = idx
                idx_to_drop = []
                idx_last = get_prev_idx(len(points_list) - 1, points_list)

                # adjust maximum distance between prev and next (if the initial distance is very high)
                total_distance = dist2prev
                dist_prev_2_next_max = min(15000, 150 + total_distance / 5)
                # print(" v1 dist_prev_2_next_max = " + str(dist_prev_2_next_max))

                # print(" idx last = " + str(idx_last))
                real_trip = False
                while dur < max_dur and idx2 < idx_last and not real_trip:
                    # current point on the path to next (which is the candidate for "back"
                    dist2next = pc.distance(points_list[idx2], points_list[get_next_idx(idx2, points_list)])
                    total_distance += dist2next

                    # update max dist (1/5 of total path distance including the candidate return)
                    dist_prev_2_next_max = min(15000, 150 + total_distance / 5)

                    # adjust maximum duration to accept later returns if the distance is very high
                    # calibrated to be :
                    # 300 at 0km
                    # 1800 at 4km
                    # max 2h30
                    max_dur = min(10000, 300 + total_distance * 3 / 8)

                    idx_to_drop.append(idx2)

                    # update next point (candidate for "back")
                    idx2 = get_next_idx(idx2, points_list)
                    point_next = points_list[idx2]

                    # update dur2base and dist between prev and next
                    dur = abs((point_next.timestamp - point.timestamp).total_seconds())
                    dist_prev_2_next = pc.distance(point_prev, point_next)

                    if not sharp_angle_real_trip(idx_to_drop, points_list):
                        # momma, we're back!! This crazy ride is over.
                        # (1) must still be within the time window
                        # (2) distance between original point and candidate should be small
                        # (3) the last leg should be large distance too (a jump) - to avoid coming back "gently"
                        comes_back = dur / max_dur + dist_prev_2_next / dist_prev_2_next_max < 2 and dist2next > dist_min
                        if debug: print(prefix + str(dur / max_dur) + " and " +
                                        str(dist_prev_2_next / dist_prev_2_next_max) + " and " + str(comes_back))
                        if comes_back:
                            break
                    else:
                        if debug: print(prefix + "abandon - real trip")
                        real_trip = True

                # drop all points between idx and idx2-1
                if comes_back:
                    idx3 = idx
                    while idx3 < idx2:
                        if debug: print(prefix + ". " + str(idx3) + ". and IDX2 = " + str(idx2))
                        points_list[idx3].drop = True
                        points_list[idx3].drop_type = 'sharp_multiple' + (go_fwd * '_fwd' + (not go_fwd) * '_not_fwd')
                        points_list[idx3].drop_only_color = debug
                        points_list[idx3].markercolor = str_color(0, 0, 255)
                        idx3 = get_next_idx(idx3, points_list)

                    # add to pts2inspect (once)
                    if debug:
                        pts_sharp.append({
                            'uidp': point.uidp,
                            'date': point.timestamp.date(),
                            'time': point.timestamp.time(),
                            'type': "sharp_mult",
                            'dist2prev': dist2prev,
                            'time2prev': time2prev,
                            'sped2prev': dist2prev / 1000 / (1 + time2prev) * 3600,
                            'oth_max_dur': max_dur,
                            'oth_max_dist': dist_prev_2_next_max,
                            'extra': points_list[idx2].timestamp.time()
                        })
                else:
                    if debug: print(prefix + "Jump that cannot be resolved + " + str(point.idx) + ", " + str(point.timestamp))

    return pts_sharp


def detect_lazy(points_list, d2p_max=150, d2n_min=750, sp2n_max=33, go_fwd=True, debug=False):
    """ detect lazy and imprecise """

    if not go_fwd:
        points_list = points_list[::-1]

    pts_lazy = []
    for idx, point in enumerate(points_list):
        if not point.drop:

            # in case we want to later drop this point (and previous ones)
            idxs_to_drop = [idx]

            # previous point, either different location or at least 2.5 minutes back
            idx_prev = get_prev_idx(idx, points_list)
            pt_prev = points_list[idx_prev]
            i = 0
            while pc.distance(point, pt_prev) < 10 and \
                  abs((point.timestamp - pt_prev.timestamp).total_seconds()) < 150 and i < 50:
                # add to the list of points to drop
                idxs_to_drop.append(idx_prev)
                idx_prev = get_prev_idx(idx_prev, points_list)
                pt_prev = points_list[idx_prev]
                i += 1

            time2prev = abs((point.timestamp - pt_prev.timestamp).total_seconds())
            dist2prev = pc.distance(point, pt_prev)

            # next point
            idx_next = get_next_idx(idx, points_list)
            pt_next = points_list[idx_next]
            time2next = abs((pt_next.timestamp - point.timestamp).total_seconds())
            dist2next = pc.distance(point, pt_next)

            # prev to next
            time_p2n = abs((pt_next.timestamp - pt_prev.timestamp).total_seconds())
            dist_p2n = pc.distance(pt_prev, pt_next)

            # lazy condition:
            # if removing current point makes the speed acceptable (from prev to next)
            if dist2prev < d2p_max \
                and dist2next > d2n_min \
                and time2next > 0 \
                and dist2next / time2next > sp2n_max > dist_p2n / time_p2n:

                    # drop ALL between idx and going back until idx_prev (not inclusive)
                    for idx_d in idxs_to_drop:
                        points_list[idx_d].drop = True
                        points_list[idx_d].drop_type = 'lazy' + (go_fwd * '_fwd' + (not go_fwd) * '_not_fwd')
                        points_list[idx_d].drop_only_color = True
                        points_list[idx_d].markercolor = str_color(0, 255, 255)

                    # add to pts2inspect - only main point
                    if debug:
                        pts_lazy.append({
                            'uidp': point.uidp,
                            'date': point.timestamp.date(),
                            'time': point.timestamp.time(),
                            'dist2next': dist2next,
                            'time2next': time2next,
                            'sped2next': dist2next / 1000 / (1 + time2next) * 3600,
                            'type': "lazy"
                        })

            # imprecise conditon
            elif dist2next > d2n_min \
                    and time2next > 0 \
                    and dist2next / time2next > sp2n_max > dist_p2n / time_p2n:

                # drop ALL between idx and going back until idx_prev (not inclusive)
                for idx_d in idxs_to_drop:
                    points_list[idx_d].drop = True
                    points_list[idx_d].drop_type = 'imprecise' + (go_fwd * '_fwd' + (not go_fwd) * '_not_fwd')
                    points_list[idx_d].drop_only_color = True
                    points_list[idx_d].markercolor = str_color(255, 255, 0)

                # add to pts2inspect - only main point
                if debug:
                    pts_lazy.append({
                        'uidp': point.uidp,
                        'date': point.timestamp.date(),
                        'time': point.timestamp.time(),
                        'dist2next': dist2next,
                        'time2next': time2next,
                        'sped2next': dist2next / 1000 / (1 + time2next) * 3600,
                        'type': "imprecise"
                    })

    return pts_lazy


def drop_first_high_speed(points_list, dist_min=750, speed_min=33.0, go_fwd=True, debug=False):
    """ Drop the first point of the day if it is followed by a high speed jump. """
    if not go_fwd:
        points_list = points_list[::-1]

    # find sequence of "first point" (within 10 m of 1st point)
    idx = 0
    if points_list[0].drop:
        idx = get_next_idx(0, points_list)
    idxs_to_drop = [idx]

    found_first_point = False
    i = 0
    while not found_first_point and i < 50:
        # if next point is close, add it to the list
        pt = points_list[idx]
        pt_next = points_list[get_next_idx(idx, points_list)]
        if pc.distance(pt, pt_next) < 10:
            idx = get_next_idx(idx, points_list)
            idxs_to_drop.append(idx)
        else:
            found_first_point = True
        i += 1

    # check if idx corresponds to a jump
    point = points_list[idx]
    # assert not point.drop

    # next point
    idx_next = get_next_idx(idx, points_list)
    pt_next = points_list[idx_next]
    time2next = abs((pt_next.timestamp - point.timestamp).total_seconds())
    dist2next = pc.distance(point, pt_next)

    # jump condition
    if (time2next == 0 and dist2next > dist_min) or \
            (time2next != 0 and dist2next > dist_min and dist2next / (1 + time2next) > speed_min):

        for idx in idxs_to_drop:
            points_list[idx].drop = True
            points_list[idx].drop_type = 'first' + (go_fwd * '_fwd' + (not go_fwd) * '_not_fwd')
            points_list[idx].drop_only_color = True
            points_list[idx].markercolor = str_color(122, 122, 122)

        # add to pts2inspect - only one point
        if debug:
            return [{
                'uidp': point.uidp,
                'date': point.timestamp.date(),
                'time': point.timestamp.time(),
                'dist2next': dist2next,
                'time2next': time2next,
                'sped2next': dist2next / 1000 / (1 + time2next) * 3600,
                'type': "high_speed",
                        }]
        else:
            return []
    return []


def detect_high_speed(points_list, dist_min=750, speed_min=33.0, go_fwd=True, debug=False):
    """ Detect high speed points. Tag only with drop_type, but not drop """

    if not go_fwd:
        points_list = points_list[::-1]

    pts_hspd = []
    for idx, point in enumerate(points_list):
        if not point.drop:
            idx_next = get_next_idx(idx, points_list)
            pt_next = points_list[idx_next]
            time2next = abs((pt_next.timestamp - point.timestamp).total_seconds())
            dist2next = pc.distance(point, pt_next)

            is_jump = False
            if time2next <= 10 and dist2next > dist_min: is_jump = True
            if time2next > 10 and \
                    above_speed_limit(dist2next, time2next, conservative=True) and \
                            dist2next > dist_min: is_jump = True

            if is_jump:
                point.drop = debug
                point.drop_type = 'unresolved' + (go_fwd * '_fwd' + (not go_fwd) * '_not_fwd')
                point.drop_only_color = True
                point.markercolor = str_color(255, 0, 0)

                # add to pts2inspect
                if debug:
                    pts_hspd.append({
                        'uidp': point.uidp,
                        'date': point.timestamp.date(),
                        'time': point.timestamp.time(),
                        'dist2next': dist2next,
                        'time2next': time2next,
                        'sped2next': dist2next / 1000 / (1 + time2next) * 3600,
                        'type': "high_speed",
                    })

    return pts_hspd


def start_location(idx, points_list, dt0=0, dt1=15, dt2=30, dt3=45, d1=100, d2=200, debug=False):
    """
    Looking backwards, start a location at idx (inclusive)
    Marks all nearby points
    """

    point = points_list[idx]

    # for all idx_1 < idx such that idx < idx_max(idx_1):
    # update score_after(idx1) upwards
    idx1 = idx - 1
    in_range = True
    while idx1 >= 0 and in_range:
        # [i_after, idx_min, idx_max] = find_points_after(idx1, points_list, dt0=0, dt1=15, dt2=30, dt3=45)
        idx_max = find_points_after(idx1, points_list, after=True, dt0=1.5, dt1=5.0, dt2=15.0)

        if idx_max >= idx:
            p = points_list[idx1]
            dist = pc.distance(point, p)
            score = (dist <= d1) + (d1 < dist < d2) * (d2 - dist) / (d2 - d1)

            # todo: change 3
            p.score_after = min(1, score + p.score_after)
            if debug:
                print([">>>", p.score_after])
        else:
            in_range = False

        idx1 -= 1


def detect_trips(points_list, debug=False):
    """
    Takes a day of data from one user, and segments it into Location/Trip/Location/Trip etc.
    """

    # assert that all points are from the same day
    dates = set(point.timestamp.date() for point in points_list)
    assert (len(dates) == 1)

    # assert the list is sorted by time
    assert (all(points_list[i] <= points_list[i + 1] for i in range(len(points_list) - 1)))

    # work on points that are not dropped
    points_list_active = [point for point in points_list if not point.drop]

    if points_list_active:
        # last point of the day
        idx = len(points_list_active) - 1
        try:
            point = points_list_active[idx]
        except Exception as e:
            print(points_list)
            raise e

        location_ongoing = False
        if point.score_before > 0.3:
            location_ongoing = True
            start_location(idx, points_list_active)
            point.trip = False

        idx = len(points_list_active) - 2
        while idx >= 0:
            try:
                point = points_list_active[idx]
            except Exception as e:
                print(points_list)
                raise e

            # for debug
            date, time = str(point.timestamp).replace("+05:30", "").split(" ")

            # if location
            if location_ongoing:

                # if still in location
                if point.score_after > 0.3:
                    point.trip = False

                # if we're starting a trip
                else:
                    if debug:
                        print("starting trip " + time)
                    point.trip = True
                    location_ongoing = False

            # if trip
            else:

                # if still trip
                if point.score_before < 0.7:
                    point.trip = True

                # if we're starting a location
                else:
                    if debug:
                        print("starting location " + time)
                    point.trip = False
                    location_ongoing = True
                    start_location(idx, points_list_active)

            idx -= 1

    # finished


