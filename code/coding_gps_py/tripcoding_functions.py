import math
import sys
import json
import os
import glob
import time
from pprint import pprint
import pointclass as pc
import utils as ut
import operator
import numpy as np
import cProfile
from joblib import Parallel, delayed
import multiprocessing
# from geopy.distance import vincenty
from geopy.distance import great_circle
from utils_geometry import diameter
from sklearn.cluster import DBSCAN
import datetime
import detect
import utils_log as ul
import logging
np.seterr(under='ignore')

"""
Example:

X orig1 =                  <- last point in previous segment
O orig  = seg [ idx ]      <-  first point in this segment  <- this will be different TYPE as orig1 (X)
O
O ...
O
O dest1 =                  < - last point in the segment
X dest  = seg [ idx + 1 ]  <- first point in NEXT segment  <- this will be different TYPE as dest1 (O)

"""

pt_central_bangalore = (12.975158, 77.590022)


def pretty_print_struct(structs, output_file):
    """ pretty print structs """
    # uidp -> date -> struct
    if structs:
        uidps = list(set([seg['uidp'] for seg in structs]))
        assert len(uidps) == 1
        uidp = uidps[0]
        struct_json = {uidp: dict()}

        current_date = ''
        for seg in structs:
            if seg['date'] != current_date:
                current_date = seg['date']
                struct_json[uidp][current_date] = []

            seg_type = "..." * (seg['trip'] == 1) + " O " * (seg['trip'] == 0)

            duration = str(math.floor(seg['dur']) % 60)
            if seg['dur'] > 60:
                duration = str(math.floor(seg['dur']/60)) + ":" + duration

            struct_json[uidp][current_date].append(
                seg_type + '| ' + seg['t_start_0'] + ' ' + seg['t_start_1'] + ' | ' + duration.rjust(5) +
                ' | ' + str(seg['sid']))
            struct_json[uidp][current_date].append('   |             |       |')
        with open(output_file, "w") as myfile:
            pprint(struct_json, myfile, width=160)

        return True
    return False


def get_prev_idx(idx, points_dict):
    # if idx not in points_dict:
    #     return -1

    keyset = {idx_new for idx_new in points_dict if idx_new < idx}

    if not keyset:
        return idx
    else:
        return max(keyset)


def get_next_idx(idx, points_dict):
    # if idx not in points_dict:
    #     return -1

    keyset = {idx_new for idx_new in points_dict if idx_new > idx}

    if not keyset:
        return idx
    else:
        return min(keyset)


def get_pointlist_seg(idx_seg, segs, point_list):

    seg = segs[idx_seg]
    i_orig = seg['start']
    i_dest = min(segs[idx_seg + 1]['start'], max([idx for idx in point_list]))

    point_list_this_seg = [point_list[idx] for idx in point_list if i_orig <= idx < i_dest]
    assert point_list_this_seg

    return point_list_this_seg


def format_time(time, seconds=False):
    """ time if float """
    if type(time) is not float:
        time = float(time)

    if seconds:
        time = int(np.round(time))
        # in seconds
        seconds = time % 60
        time = int((time - seconds) / 60)
        minutes = time % 60
        hours = int((time - minutes) / 60)
        return str(hours) + ":" + str(minutes).zfill(2) + ":" + str(seconds).zfill(2)
    else:
        hour = int(np.floor(time))
        minutes = int(np.round((time - hour) * 60))
        return str(hour) + ":" + str(minutes).zfill(2)


def minute_hhmm(tt, seconds=False):
    # compute minutes from midnight given 5 digit string hh:mm
    assert tt[2] == ":"
    hh = int(tt[0:2])
    mm = int(tt[3:5])
    if seconds:
        assert tt[2] == tt[5] == ":"
        ss = int(tt[6:8])
        return 3600 * hh + 60 * mm + ss
    else:
        return 60 * hh + mm


def minutes_hhmmss(tt):
    # compute seconds from midnight given 8 digit hh:mm:ss
    assert tt[2] == tt[5] == ":"
    hh = int(tt[0:2])
    mm = int(tt[3:5])


def dur_hhmm(orig, dest):
    orig_hh = int(orig[0:2])
    orig_mm = int(orig[3:5])
    dest_hh = int(dest[0:2])
    dest_mm = int(dest[3:5])

    return abs(float(dest_mm - orig_mm) + float(dest_hh - orig_hh) * 60)


def dur_hhmm_points(orig_time, dest_time):
    orig_hh = int(orig_time[0:2])
    orig_mm = int(orig_time[3:5])
    orig_ss = int(orig_time[6:8])
    # fractional minutes
    orig_mf = float(orig_mm) + float(orig_hh) * 60 + float(orig_ss) / 60
    orig_mf = round(orig_mf, 2)

    dest_hh = int(dest_time[0:2])
    dest_mm = int(dest_time[3:5])
    dest_ss = int(dest_time[6:8])
    # fractional minutes
    dest_mf = float(dest_mm) + float(dest_hh) * 60 + float(dest_ss) / 60
    dest_mf = round(dest_mf, 2)

    return abs(dest_mf - orig_mf)


def segs_stats_ids(segs):
    """ add ids """
    for idx_seg, seg in enumerate(segs):
        seg['id'] = idx_seg
    return True


def segs_stats_cluster_codes(segs):
    # update trip and jump origin and destinations
    # go forward - origins
    loc_cluster = -1
    loc_lat = -1
    loc_lon = -1
    jump_ongoing = False
    for seg in segs:
        if seg['type'] == 'loc':
            loc_cluster = seg.get('c_orig', '')
            loc_lat = seg.get('center_lat', '')
            loc_lon = seg.get('center_lon', '')
            jump_ongoing = False
        elif seg['type'] == 'trip':
            seg['c_orig'] = loc_cluster
            seg['c_orig_lat'] = loc_lat
            seg['c_orig_lon'] = loc_lon
            loc_cluster = -1
            loc_lat = -1
            loc_lon = -1
            if jump_ongoing:
                seg['jump_pre'] = 1
            else:
                seg['jump_pre'] = 0
            jump_ongoing = False
        elif seg['type'] == 'jump' or seg['type'] == 'jumptrip':
            if seg.get('duration', 31) > 30 or seg.get('line_dist', 5001) > 5000:
                jump_ongoing = True
        elif seg['type'] == 'gap':
            jump_ongoing = False
            seg['c_orig'] = loc_cluster
            seg['c_orig_lat'] = loc_lat
            seg['c_orig_lon'] = loc_lon
        elif seg['type'] == 'end':
            pass
        else:
            print("unexpected seg type - forward")
            raise ArithmeticError

    # go backwards - destinations
    loc_cluster = -1
    loc_lat = -1
    loc_lon = -1
    jump_ongoing = False
    for seg in reversed(segs):
        if seg['type'] == 'loc':
            loc_cluster = seg.get('c_orig', '')
            loc_lat = seg.get('center_lat', '')
            loc_lon = seg.get('center_lon', '')
            jump_ongoing = False
        elif seg['type'] == 'trip':
            seg['c_dest'] = loc_cluster
            seg['c_dest_lat'] = loc_lat
            seg['c_dest_lon'] = loc_lon
            loc_cluster = -1
            loc_lat = -1
            loc_lon = -1
            if jump_ongoing:
                seg['jump_post'] = 1
            else:
                seg['jump_post'] = 0
            jump_ongoing = False
        elif seg['type'] == 'jump' or seg['type'] == 'jumptrip':
            if seg.get('duration', 31) > 30 or seg.get('line_dist', 5001) > 5000:
                jump_ongoing = True
        elif seg['type'] == 'gap':
            jump_ongoing = False
            seg['c_orig'] = loc_cluster
            seg['c_orig_lat'] = loc_lat
            seg['c_orig_lon'] = loc_lon
        elif seg['type'] == 'end':
            pass
        else:
            print("unexpected seg type - backward")
            raise ArithmeticError
    return True


def segs_stats(segs, point_dict):

    for idx_seg, seg in enumerate(segs[:-1]):
        try:
            max_idx = max(point_dict)

            i_orig = seg['start']
            i_orig1 = get_prev_idx(i_orig, point_dict)
            i_dest = min(segs[idx_seg + 1]['start'], max_idx)
            i_dest1 = get_prev_idx(i_dest, point_dict)

            assert i_dest > i_orig or (i_dest == max_idx)

            orig = point_dict[i_orig]
            dest = point_dict[i_dest]

            # arrival and departure times
            seg['t00'] = point_dict[i_orig1]['time'][0:5]
            seg['t01'] = orig['time'][0:5]
            seg['t10'] = point_dict[i_dest1]['time'][0:5]
            seg['t11'] = dest['time'][0:5]

            seg['t00_'] = minute_hhmm(point_dict[i_orig1]['time'], seconds=True)
            seg['t01_'] = minute_hhmm(orig['time'], seconds=True)
            seg['t10_'] = minute_hhmm(point_dict[i_dest1]['time'], seconds=True)
            seg['t11_'] = minute_hhmm(dest['time'], seconds=True)

            # best guess arrival and departure times
            if seg['type'] == 'loc':
                # start location exactly when it starts
                seg['t0bg_'] = seg['t01_']
                seg['t0bg'] = format_time(seg['t0bg_'], seconds=True)
            elif seg['type'] == 'trip':
                # start trip earlier but not more than 5 mins from start
                seg['t0bg_'] = max(seg['t00_'], seg['t01_'] - 5 * 60)
                seg['t0bg'] = format_time(seg['t0bg_'], seconds=True)
            else:
                # avearge of jump or gap but not too far from the origin (5 minutes)
                # seg['t0bg_'] = min((minute_hhmm(seg['t01']) + minute_hhmm(seg['t11'])) * 0.5,
                #                    minute_hhmm(seg['t01']) + 5)
                # exactly origin
                seg['t0bg_'] = seg['t01_']
                seg['t0bg'] = format_time(seg['t0bg_'], seconds=True)


            # duration
            seg['duration'] = dur_hhmm_points(orig['time'], dest['time'])
            seg['line_dist'] = pc.distance_dict(orig, dest)

            # number of points (including endpoint)
            seg['np'] = len([idx for idx in point_dict if i_dest >= idx >= i_orig])
            assert seg['np'] >= 1

            # define list of points associated with this segment
            # all points between i_orig and i_dest
            # point_list_this_seg = [point_dict[idx] for idx in point_dict if i_orig <= idx < max(i_orig + 1, i_dest)]
            # # if seg['type'] == 'jumptrip':
            # if seg['type'] in ['jumptrip', 'jump', 'trip']:
            #     point_list_this_seg = [point_dict[idx] for idx in point_dict if i_orig <= idx <= i_dest]
            point_list_this_seg = [point_dict[idx] for idx in sorted(point_dict.keys()) if i_orig <= idx <= i_dest]
            point_list_this_seg_proper = point_list_this_seg[:-1] or point_list_this_seg
            # print(seg['t01'])
            # print(i_orig)
            # print(i_dest)
            # pprint(point_list_this_seg)

            # mark distance to central Bangalore
            # pt_central_bangalore = (12.975158, 77.590022)
            for idx, pt in enumerate(point_list_this_seg):
                pt_latlong = (float(pt['latitude']), float(pt['longitude']))
                pt['dist2c'] = great_circle(pt_latlong, pt_central_bangalore).meters
                pt['outstation'] = pt['dist2c'] > 18000
                # todo : include here radius to home location (if out of station)

            # define location centre
            if seg['type'] == 'loc':
                seg['center_lat'], seg['center_lon'], seg['center_acc'], seg['center_std'] = \
                        pc.point_list_center(point_list_this_seg_proper, accuracy_fraction=80, min_acc=0)

            # default outstation
            seg['outstation'] = 0
            seg['outstation_partial'] = 0
            seg['outstation_share'] = 0.0
            seg['outstation_full'] = 0
            seg['outstation_np'] = 0

            # outstation vars
            if seg['np'] > 0:
                # vector of distances
                dist_list = np.array([pt['dist2c'] for pt in point_list_this_seg_proper])

                # outstation
                if point_list_this_seg_proper:
                    seg['d2b_min'] = np.min(dist_list)
                    seg['d2b_p20'] = np.percentile(dist_list, 20)
                    seg['d2b_avg'] = np.average(dist_list)
                    seg['d2b_p80'] = np.percentile(dist_list, 80)
                    seg['d2b_max'] = np.max(dist_list)
                    seg['outstation_np'] = np.sum(dist_list > 18000)

                    # make sure that all jumps have exactly two points
                    # if seg['type'] in ['jumptrip', 'jump']:
                    #     if seg['np'] != 2:
                    #         pprint(seg)
                    #     assert seg['np'] == 2

                    if seg['d2b_p20'] > 18000 or (seg['type'] in ['jumptrip', 'jump'] and seg['d2b_max'] >= 18000):
                        seg['outstation'] = 1
                        seg['outstation_full'] = 1
                    elif seg['d2b_p20'] <= 18000 < seg['d2b_p80']:
                        seg['outstation'] = 1
                        seg['outstation_partial'] = 1
                    else:  # seg['d2b_p80'] <= 18000
                        pass

            """ distance in bangalore """
            # define path length, and path length inside Bangalore
            if seg['type'] in ['trip', 'jump', 'jumptrip'] and not seg.get('subtrips_ids', False):
                seg['path_length'] = sum([pt['dist2next'] for pt in point_list_this_seg_proper])

                # add first bit of the trip
                if seg['type'] == 'trip':
                    seg['path_length'] += point_list_this_seg_proper[0]['dist2prev']

                # sum up distance over all pairs of consecutive points, BOTH NOT outstation
                seg['path_length_inbang'] = sum([point_list_this_seg[i]['dist2next']
                                                 for i in range(len(point_list_this_seg) - 1)
                                                 if not point_list_this_seg[i]['outstation']
                                                 and not point_list_this_seg[i + 1]['outstation']])

                # add first bit of the trip (if both endpoints are inside Bangalore)
                idx0 = int(point_list_this_seg[0]['idx'])
                idxm1 = get_prev_idx(idx0, point_dict)
                if seg['type'] == 'trip' and not point_dict[idx0]['outstation'] and not point_dict[idxm1]['outstation']:
                    seg['path_length_inbang'] += point_dict[idx0]['dist2prev']
            # done

        except Exception as e:
            print("DEBUG:")
            # pprint(segs)
            pprint(seg)
            # pprint(segs[idx_seg + 1])
            print(i_orig)
            print(len(point_dict))
            raise e

    # update best guess estimate of trip distance and duration
    for idx_seg, seg in enumerate(segs[:-1]):
        if seg['type'] == 'trip':
            seg['path_length_inbang_tot'] = seg['path_length_inbang']
            seg['path_length_tot'] = seg['path_length']
            seg['duration_tot'] = seg['duration']

            # todo: but only if not outstation jump
            if idx_seg > 0:
                if segs[idx_seg - 1]['type'] == 'jump' and not segs[idx_seg - 1]['outstation_full']:
                    seg['path_length_inbang_tot'] += segs[idx_seg - 1]['line_dist']
                    seg['path_length_tot'] += segs[idx_seg - 1]['line_dist']
                    seg['duration_tot'] += segs[idx_seg - 1]['duration']

            # todo: but only if not outstation jump
            if segs[idx_seg + 1]['type'] == 'jumptrip' and not segs[idx_seg + 1]['outstation_full']:
                seg['path_length_inbang_tot'] += segs[idx_seg + 1]['line_dist']
                seg['path_length_tot'] += segs[idx_seg + 1]['line_dist']
                seg['duration_tot'] += segs[idx_seg + 1]['duration']

            if seg['duration_tot'] > 0:
                seg['speed_kmph'] = (seg['path_length_tot'] / 1000) / (seg['duration_tot'] / 60)
            else:
                seg['speed_kmph'] = -1

    return True


def point_simple_stats(points_list, point_dict, distances_removed=None):

    # times
    for point in points_list:
        point['hh'] = int(point['time'][0:2])
        point['mm'] = int(point['time'][3:5])
        point['ss'] = int(point['time'][6:8])

        # fractional minutes
        point['mf'] = float(point['mm']) + float(point['hh']) * 60 + float(point['ss']) / 60
        point['mf'] = round(point['mf'], 2)

    for idx, point in enumerate(points_list):
        pt_prev = point
        if idx > 0:
            pt_prev = points_list[idx - 1]
        pt_next = point
        if idx + 1 < len(points_list):
            pt_next = points_list[idx + 1]

        # times
        point['time2next'] = dur_hhmm_points(point['time'], pt_next['time'])
        point['time2prev'] = dur_hhmm_points(point['time'], pt_prev['time'])

        point['dist2prev'] = great_circle((float(point['latitude']), float(point['longitude'])),
                                          (float(pt_prev['latitude']), float(pt_prev['longitude']))).meters
        point['dist2next'] = great_circle((float(point['latitude']), float(point['longitude'])),
                                          (float(pt_next['latitude']), float(pt_next['longitude']))).meters

        """ filter out unreasonable speed """
        # 750 meters minimum and above (conservative) speed limit
        for suf in ['next', 'prev']:
            dist = point['dist2' + suf]  # meters
            dur = point['time2' + suf] * 60  # seconds

            if (dur == 0 and dist > 750) or \
               (dur != 0 and dist > 750 and detect.above_speed_limit(dist=dist, dur=dur, conservative=True)):

                # add to list
                if distances_removed is not None:
                    # distance to central Bangalore
                    dist2bang = great_circle((float(point['latitude']), float(point['longitude'])),
                                             pt_central_bangalore).meters
                    distances_removed.append({
                        'uidp': point['uidp'],
                        'date': point['date'],
                        'time': point['time'],
                        'suf': suf,
                        'time2': point['time2' + suf],
                        'dist2': point['dist2' + suf],
                        'dist2bang': dist2bang
                    })
                # point['dist2' + suf] = 0 # todo: DO NOT DROP for now

        # refresh in dictionary
        if point_dict:
            point_dict[int(point['idx'])]['dist2prev'] = point['dist2prev']
            point_dict[int(point['idx'])]['dist2next'] = point['dist2next']
            point_dict[int(point['idx'])]['time2prev'] = point['time2prev']
            point_dict[int(point['idx'])]['time2next'] = point['time2next']

    # return points_list
    return True


def compute_times_struct(idx_seg, seg_struct, segs_dict, segs_struct, point_dict):
    """
    START
        for normal segments: start_0 = t00 (=last point before segment)
        for normal segments: start_1 = t01 (=first point on the trip/ at location)
        for jumps/gaps/etc, t01==t10 is the start of the jump, and t11 is the end, so use these.

    To pick the *best* estimate:
        't0bg_'
    STOP
        always t10 and t11
    """
    """ distances """
    # start and stop distance covered between upper and lower bound
    if seg_struct['type_first'] in ['trip', 'loc']:
        seg_struct['dist_start'] = point_dict[segs_dict[seg_struct['segs'][0]]['start']]['dist2prev']
    else:
        # for gap, jump and jumptrip
        seg_struct['dist_start'] = point_dict[segs_dict[seg_struct['segs'][0]]['start']]['dist2next']
    # it's simpler at the end of the trip
    seg_struct['dist_stop'] = point_dict[segs_dict[seg_struct['segs'][-1]]['start']]['dist2next']

    """ times """
    # upper and lower bound for departure (start) time - in SECONDS
    seg_struct['t_start_0_'] = seg_struct['t00_']
    seg_struct['t_start_1_'] = seg_struct['t01_']
    # different for jumps and gaps - take start and finish - in SECONDS
    if seg_struct['type_first'] in ['gap', 'jump', 'jumptrip']:
        seg_struct['t_start_0_'] = seg_struct['t01_']
        seg_struct['t_start_1_'] = segs_dict[seg_struct['segs'][0]]['t11_']

    # ... and for arrival - in SECONDS
    seg_struct['t_stop_0_'] = segs_dict[seg_struct['segs'][-1]]['t10_']
    seg_struct['t_stop_1_'] = segs_dict[seg_struct['segs'][-1]]['t11_']

    # format as nice time
    for field in ['t_start_0', 't_start_1', 't_stop_0', 't_stop_1']:
        seg_struct[field] = format_time(seg_struct[field + '_'], seconds=True)

    # best guess:
    # seg_struct['t_start_bg'] = seg_struct['t0bg_']
    if idx_seg < len(segs_struct) - 1:
        seg_struct['t_stop_bg'] = segs_struct[idx_seg + 1]['t_start_bg']
        seg_struct['t_stop_bg_'] = segs_struct[idx_seg + 1]['t_start_bg_']
    else:
        seg_struct['t_stop_bg'] = "23:59:59"
        seg_struct['t_stop_bg_'] = 24 * 3600 - 1
    seg_struct['dur_bg_'] = seg_struct['t_stop_bg_'] - seg_struct['t_start_bg_']
    seg_struct['dur_bg_mm'] = np.floor(seg_struct['dur_bg_'] / 6) / 10

    # departure/arrival precision - in SECONDS
    seg_struct['t_start_prec'] = seg_struct['t_start_1_'] - seg_struct['t_start_0_']
    seg_struct['t_stop_prec'] = seg_struct['t_stop_1_'] - seg_struct['t_stop_0_']

    """ durations in minutes """
    seg_struct['dur'] = sum([segs_dict[idx]['duration'] for idx in seg_struct['segs']])

    seg_struct['dur_trip'] = sum([segs_dict[idx]['duration'] for idx in seg_struct['segs']
                                  if segs_dict[idx]['type'] in ['trip']])

    seg_struct['dur_jump'] = sum([segs_dict[idx]['duration'] for idx in seg_struct['segs']
                                  if segs_dict[idx]['type'] in ['jump', 'jumptrip']])

    seg_struct['dur_loc'] = sum([segs_dict[idx]['duration'] for idx in seg_struct['segs']
                                 if segs_dict[idx]['type'] == 'loc'])

    seg_struct['dur_gap'] = sum([segs_dict[idx]['duration'] for idx in seg_struct['segs']
                                 if segs_dict[idx]['type'] == 'gap'])


def print_segs(segs, point_list):
    for seg in segs[:-1]:
        try:
            pprint(seg)
            # print("start:" + point_list[seg['start']]['time'] + ": " + seg['type'] + ": idx:" + str(seg['start']))
        except Exception as e:
            print(seg)
            print(segs)
            raise e


def find_segs_naive(points_list):
    """
    a trip is all the points with trip==1
    trips and locations are separated by "bridge" links
    :param points_list:
    :return:
    """

    segs = []

    my_date = points_list[0]['date']

    # first points
    assert points_list[0]['trip'] in ['0', '1']
    trip = int(points_list[0]['trip'])

    current_seg = {
        'date': my_date,
        'start': min([int(point['idx']) for point in points_list]),
        'type': 'trip' * trip + 'loc' * (1 - trip)
    }

    # go through and accumulate loc/trip
    for idx, point in enumerate(points_list):
        assert point.get('drop', '0') != '1'
        assert point['trip'] in ['0', '1']
        trip = int(point['trip'])
        # ongoing trip
        if current_seg['type'] == 'trip':
            # continue
            if trip:
                pass
            # end current trip and start new location
            else:
                segs.append(current_seg)
                current_seg = {
                    'date': my_date,
                    'start': int(point['idx']),
                    'type': 'loc'
                }
        # ongoing location
        else:
            # continue
            if not trip:
                pass
            # end current location and start new trip
            else:
                segs.append(current_seg)
                current_seg = {
                    'date': my_date,
                    'start': int(point['idx']),
                    'type': 'trip'
                }

    segs.append(current_seg)

    # add token end segment
    current_seg = {
        'date': my_date,
        'start': int(points_list[-1]['idx']) + 1,
        'type': 'end'
    }
    segs.append(current_seg)

    return segs


def process_ok_low_acc(segs, point_dict, points_dict_full, logger=None, silent=False):
    for idx_seg, seg in enumerate(segs):
        if seg['type'] == 'loc':
            # idx_start = seg['start']
            # idx_end_full = get_prev_idx(segs[idx_seg + 1]['start'], points_dict_full)

            # idx_start = first full point after end of previous location
            idx_start = get_prev_idx(seg['start'], point_dict)
            idx_start = min(get_next_idx(idx_start, points_dict_full), seg['start'])

            idx_end_full = segs[idx_seg + 1]['start']

            # not for degenerate locations
            # if idx_start < idx_end_full:

            # get location centre coordinates
            idx_list = [idx for idx in range(idx_start, idx_end_full) if idx in point_dict]
            point_list_this_seg = [point_dict[idx] for idx in idx_list]
            try:
                center_lat, center_lon, _, _ = \
                pc.point_list_center(point_list_this_seg, accuracy_fraction=80, min_acc=0)
            except Exception as e:
                pprint(segs)
                print(seg)
                print(idx_start)
                print(idx_end_full)
                raise e

            # include low accuracy points that are within range
            idx_list_full = [idx for idx in range(idx_start, idx_end_full) if idx in points_dict_full]
            idx_list_full = [idx for idx in idx_list_full if
                             great_circle((center_lat, center_lon),
                                          (float(points_dict_full[idx]['latitude']),
                                       float(points_dict_full[idx]['longitude']))).meters
                             < min(200, float(points_dict_full[idx]['accuracy_level']))
                             and (points_dict_full[idx]['drop'] == '0'
                                  or points_dict_full[idx].get('drop_type', '') == 'low_acc')]

            # filter: at most one point per minute
            idx_list_temp = idx_list_full
            idx_list_full = []
            for i, idx in enumerate(idx_list_temp):
                if i == 0:
                    idx_list_full.append(idx)
                else:
                    point = points_dict_full[idx]
                    idx_last = idx_list_temp[i - 1]
                    pt_last = points_dict_full[idx_last]
                    time2last = dur_hhmm_points(point['time'], pt_last['time'])
                    if time2last >= 1.0:
                        idx_list_full.append(idx)

            # change start point if necessary
            start_new = min(set(idx_list_full).union(set(idx_list)))
            assert start_new <= seg['start']
            seg['start'] = start_new

            # add to list of points
            for idx in idx_list_full:
                if idx not in point_dict:
                    if not silent:
                        logger.warning("Added low acc point " +
                                   points_dict_full[idx]['date'] + " " + points_dict_full[idx]['time'] +
                                   " acc=" + points_dict_full[idx]['accuracy_level'])
                    point_dict[idx] = points_dict_full[idx]
                    point_dict[idx]['dist2next'] = -1.0
                    point_dict[idx]['dist2prev'] = -1.0
                    assert point_dict[idx]['drop'] == '1'


def process_seg_gap(idx_seg, segs_new, points_list, logger=None, silent=False):
    """ insert a single gap/jump/jumptrip at a certain index """
    idx0 = min(points_list.keys())
    my_date = points_list[idx0]['date']

    seg = segs_new[idx_seg]
    # check if there is any jump or gap
    if seg['type'] == 'loc':
        # LOC - inside
        # idx = seg['start']
        # idx_next = get_prev_idx(segs_new[idx_seg + 1]['start'], points_list)
        # while idx < idx_next:

        idx_start = seg['start']
        idx_end = get_prev_idx(segs_new[idx_seg + 1]['start'], points_list)
        idx_list = [idx for idx in range(idx_start, idx_end) if idx in points_list]
        for idx in idx_list:
            point = points_list[idx]
            pt_next = points_list[get_next_idx(idx, points_list)]

            # LOC - if time to next is large
            delta_t = dur_hhmm_points(point['time'], pt_next['time'])
            assert math.isclose(delta_t, point['time2next'])
            if delta_t > 60:
                # LOC - insert gap
                if not silent: logger.warning("SEG GAP (delta_t > 60) LOC - insert gap after segment " + str(idx_seg))
                newseg_gap = {
                    'date': my_date,
                    'type': 'gap',
                    'start': idx
                }
                newseg_nextloc = {
                    'date': my_date,
                    'type': 'loc',
                    'start': get_next_idx(idx, points_list)
                }

                # gap (replace old loc) - loc
                if idx == idx_list[0]:
                    assert idx_list[0] == seg['start']
                    segs_new[idx_seg] = newseg_gap
                    segs_new.insert(idx_seg + 1, newseg_nextloc)
                    idx_seg += 1
                    return idx_seg
                # loc - gap (new) - loc (new)
                else:
                    segs_new.insert(idx_seg + 1, newseg_gap)
                    segs_new.insert(idx_seg + 2, newseg_nextloc)
                    idx_seg += 2
                    return idx_seg

        # LOC - last point (connection with next segment)
        idx = idx_end
        assert idx >= seg['start']
        point = points_list[idx]
        delta_t = point['time2next']
        delta_d = point['dist2next']

        if delta_d > 1000 or (delta_d > 500 and delta_t > 10):
            newseg_jump = {
                'date': my_date,
                'type': 'jump',
                'start': idx
            }

            if idx == seg['start']:
                # LOC - jump - replace exising location
                if not silent: logger.warning("SEG GAP (delta_d > 1000 or (delta_d > 500 and delta_t > 10) LOC - jump - replace exising location) after segment " + str(idx_seg))
                segs_new[idx_seg] = newseg_jump
                idx_seg += 1

                # pprint(segs_new[:(idx_seg + 1)])
                # print("reached this point... jump replace")
                # input()

                return idx_seg
            else:
                # LOC - jump
                if not silent: logger.warning("SEG GAP (delta_d > 1000 or (delta_d > 500 and delta_t > 10) LOC - jump) after segment " + str(idx_seg))
                segs_new.insert(idx_seg + 1, newseg_jump)
                idx_seg += 2

                return idx_seg

        elif delta_t > 20:
            newseg_gap = {
                'date': my_date,
                'type': 'gap',
                'start': idx
            }
            if idx == seg['start']:
                # LOC - gap - replace location with this
                if not silent: logger.warning("SEG GAP (delta_t > 20) LOC - gap - replace location with this) after segment " + str(idx_seg))
                segs_new[idx_seg] = newseg_gap
                idx_seg += 1

                # pprint(segs_new[:(idx_seg + 1)])
                # print("reached this point... gap replace ")
                # input()

                return idx_seg
            else:
                # LOC - gap
                if not silent: logger.warning("SEG GAP (delta_t > 20) LOC - gap) after segment " + str(idx_seg))
                segs_new.insert(idx_seg + 1, newseg_gap)
                idx_seg += 2
                return idx_seg

        # no new segments
        idx_seg += 1
        return idx_seg

    else:
        assert seg['type'] == 'trip'
        # TRIP - all points
        # idx = seg['start']
        # idx_next = segs_new[idx_seg + 1]['start']
        # while idx < idx_next:

        idx_list = [idx for idx in range(seg['start'], segs_new[idx_seg + 1]['start']) if idx in points_list]
        idx_end = idx_list[-1]
        for idx in idx_list:
            point = points_list[idx]
            delta_t = point['time2next']
            delta_d = point['dist2next']

            # TRIP - jumptrip if (v large distance) or large time and medium distance
            if delta_d > 2000 or \
                    (delta_d > 500 and delta_t > 15) or \
                    (delta_d > 300 and delta_t > 20):
                # TRIP - insert trip_jump
                newseg_jump = {
                    'date': my_date,
                    'type': 'jumptrip',
                    'start': idx
                }
                newseg_nexttrip = {
                    'date': my_date,
                    'type': 'trip',
                    'start': get_next_idx(idx, points_list)
                }

                if seg['start'] == idx == idx_end:
                    """ replace trip with jumptrip """
                    if not silent:
                        logger.warning("SEG GAP (d > 2000 or (d > 500 and t > 15) or (d > 300 and t > 20)) "
                                   "TRIP replace trip with jumptrip. segment " + str(idx_seg))
                    segs_new[idx_seg] = newseg_jump

                    # pprint(segs_new[:(idx_seg + 1)])
                    # print("reached this point... case A")
                    # input()

                    idx_seg += 1
                    return idx_seg

                elif seg['start'] == idx < idx_end:
                    """ replace trip with jumptrip-trip """
                    if not silent:
                        logger.warning("SEG GAP (d > 2000 or (d > 500 and t > 15) or (d > 300 and t > 20)) "
                                   "TRIP replace trip with jumptrip-trip. segment " + str(idx_seg))
                    segs_new[idx_seg] = newseg_jump
                    segs_new.insert(idx_seg + 1, newseg_nexttrip)

                    # pprint(segs_new[:(idx_seg + 2)])
                    # print("reached this point... case B")
                    # input()

                    idx_seg += 1
                    return idx_seg

                elif seg['start'] < idx == idx_end:
                    """ replace trip with trip-jumptrip """
                    if not silent:
                        logger.warning("SEG GAP (d > 2000 or (d > 500 and t > 15) or (d > 300 and t > 20)) "
                                   "TRIP replace trip with trip-jumptrip. segment " + str(idx_seg))
                    segs_new.insert(idx_seg + 1, newseg_jump)

                    # pprint(segs_new[:(idx_seg + 2)])
                    # print("reached this point... case C")
                    # input()

                    idx_seg += 2
                    return idx_seg
                else:
                    """ replace trip with trip-jumptrip-trip """
                    if not silent:
                        logger.warning("SEG GAP (d > 2000 or (d > 500 and t > 15) or (d > 300 and t > 20)) "
                                   "TRIP replace trip with trip-jumptrip-trip. segment " + str(idx_seg))
                    assert seg['start'] < idx < idx_end
                    segs_new.insert(idx_seg + 1, newseg_jump)
                    segs_new.insert(idx_seg + 2, newseg_nexttrip)

                    # pprint(segs_new[:(idx_seg + 2)])
                    # print("reached this point... case D")
                    # input()

                    idx_seg += 2
                    return idx_seg

            # INSERT LOCATION if the distance is small and there is a time gap
            elif delta_t > 10 and delta_d < 300:

                newseg_loc = {
                    'date': my_date,
                    'type': 'loc',
                    'start': idx
                }
                newseg_gap = {
                    'date': my_date,
                    'type': 'gap',
                    'start': idx
                }
                newseg_nexttrip = {
                    'date': my_date,
                    'type': 'trip',
                    'start': get_next_idx(idx, points_list)
                }

                if seg['start'] == idx == idx_end:
                    """ replace trip with gap """
                    if not silent:
                        logger.warning("SEG GAP (t > 10 and d < 300) TRIP replace trip with gap. segment " + str(idx_seg))
                    segs_new[idx_seg] = newseg_gap

                    # pprint(segs_new[:(idx_seg + 1)])
                    # print("reached this point... case A")
                    # input()

                    idx_seg += 1
                    return idx_seg

                elif seg['start'] == idx < idx_end:
                    """ replace trip with loc-trip """
                    if not silent:
                        logger.warning("SEG GAP (t > 10 and d < 300) TRIP replace trip with loc-trip. segment " + str(idx_seg))
                    segs_new[idx_seg] = newseg_loc
                    segs_new.insert(idx_seg + 1, newseg_nexttrip)

                    # pprint(segs_new[:(idx_seg + 2)])
                    # print("reached this point... case B")
                    # input()

                    idx_seg += 1
                    return idx_seg

                elif seg['start'] < idx == idx_end:
                    """ replace trip with trip-gap """
                    if not silent:
                        logger.warning("SEG GAP (t > 10 and d < 300) TRIP replace trip with trip-gap. segment " + str(idx_seg))
                    segs_new.insert(idx_seg + 1, newseg_gap)

                    # pprint(segs_new[:(idx_seg + 2)])
                    # print("reached this point... case C")
                    # input()

                    idx_seg += 2
                    return idx_seg
                else:
                    """ replace trip with trip-loc-trip """
                    if not silent:
                        logger.warning("SEG GAP (t > 10 and d < 300) TRIP replace trip with trip-loc-trip. segment " + str(idx_seg))
                    assert seg['start'] < idx < idx_end
                    segs_new.insert(idx_seg + 1, newseg_loc)
                    segs_new.insert(idx_seg + 2, newseg_nexttrip)

                    # pprint(segs_new[:(idx_seg + 2)])
                    # print("reached this point... case D")
                    # input()

                    idx_seg += 2
                    return idx_seg

        # no new segments
        idx_seg += 1
        return idx_seg


def process_segs_gaps(segs, points_list, logger=None, silent=False):
    """ go through list of segments and insert gaps, jumps and jumptrips where needed """

    segs_new = segs
    idx_seg = 0
    while idx_seg < len(segs_new) - 1:
        if not silent: print('.', end='', flush=True)
        idx_seg = process_seg_gap(idx_seg, segs_new, points_list, logger=logger, silent=silent)

    if not silent:
        print('.')
        logger.info('.')
    # add the end timestamp and exit
    return segs_new


def process_short_trip(idx_seg, segs_new, points_list, logger=None, silent=False):

    # previous segment exists
    assert idx_seg > 0
    # there is a next segment
    assert segs_new[-1]['type'] == 'end'

    if idx_seg + 1 < len(segs_new) - 1:
        seg_m1 = segs_new[idx_seg - 1]
        seg = segs_new[idx_seg]
        seg_p1 = segs_new[idx_seg + 1]

        # short trip/jump sandwiched between locations
        if seg_m1['type'] == seg_p1['type'] == 'loc' and seg['type'] in ['trip', 'jump', 'jumptrip']:
            duration = seg['duration']
            pm1 = (seg_m1['center_lat'], seg_m1['center_lon'])
            pp1 = (seg_p1['center_lat'], seg_p1['center_lon'])
            distance = great_circle(pm1, pp1).meters

            # compute trip diameter (furthest away points)
            trip_point_list = get_pointlist_seg(idx_seg, segs_new, points_list)
            trip_point_list = [(float(point['latitude']), float(point['longitude'])) for point in trip_point_list]
            if len(trip_point_list) > 1:
                trip_diameter, _ = diameter(trip_point_list)
            else:
                # assume diameter when we only have 1 point
                assert len(trip_point_list) == 1
                trip_diameter = seg['line_dist']
            # print("Diameter " + str(trip_diameter))

            # todo: make this condition much more stringent - only if the locations are the same,
            # todo: distance and time are small.
            if ((distance < 400 and duration < 10) or (distance < 200 and duration < 20)) and trip_diameter < 300:
                if not silent: logger.warning("REMOVE SHORT TRIP " + str(idx_seg) + " " + seg['t01'])
                # pop trip and merge the two locations
                del segs_new[idx_seg]
                del segs_new[idx_seg]
                # update all stats
                # segs_stats(segs_new, points_list)
                # process the same (used to be idx_seg + 2
                return idx_seg

        # mini trip
        elif seg['type'] == 'trip' and seg['duration'] < 1 and seg['line_dist'] < 500:
            if not silent: logger.warning("REMOVE SHORT TRIP < 1 minute" + str(idx_seg))
            del segs_new[idx_seg]
            # update all stats
            # segs_stats(segs_new, points_list)
            # process the same (used to be idx_seg + 2
            return idx_seg

    return idx_seg + 1


def process_short_trips(segs, points_list, logger=None, silent=False):

    # looking for short trips
    segs_new = segs
    idx_seg = 1
    while idx_seg < len(segs_new) - 1:
        idx_seg = process_short_trip(idx_seg, segs_new, points_list, logger=logger, silent=silent)

    # refresh
    segs_stats(segs_new, points_list)
    return segs_new


def find_next(idx_seg, segs, seg_type):
    idx_seg += 1
    while idx_seg < len(segs):
        if segs[idx_seg]['type'] == seg_type:
            return idx_seg
        else:
            idx_seg += 1
    return 0


def process_loc_clusters(uidp, segs, PATHOUT, logger=None, silent=False):

    locs = [seg for seg in segs if seg['type'] == 'loc']
    lats = [float(loc['center_lat']) for loc in locs]
    lons = [float(loc['center_lon']) for loc in locs]

    # coords = np.matrix([lats, lons]).transpose()
    coords = np.asarray([lats, lons]).transpose()

    kms_per_radian = 6371.0088
    epsilon = 0.3 / kms_per_radian
    db = DBSCAN(eps=epsilon, min_samples=1, algorithm='ball_tree', metric='haversine').fit(np.radians(coords))

    num_clusters = len(set(db.labels_))
    if not silent:
        logger.info("number of locations: " + str(len(coords)))
        logger.info("number of clusters: " + str(num_clusters))

    # assign cluster labels to segments
    for idx, seg in zip(list(db.labels_), locs):
        seg['c_orig'] = idx

    # collect clusters
    clusters = [{'point_list': [], 'id_old': i} for i in range(num_clusters)]
    for idx, loc in enumerate(locs):
        idx_cluster = loc['c_orig']
        clusters[idx_cluster]['point_list'].append({
            'latitude': lats[idx],
            'longitude': lons[idx],
            'np': loc['np']
        })

    # compute clusters stats
    for cluster in clusters:
        cluster['np'] = sum([pt['np'] for pt in cluster['point_list']])
        cluster['nloc'] = len(cluster['point_list'])
        cluster['latitude'], cluster['longitude'], cluster['accuracy_level'], cluster['center_std'] = \
            pc.point_list_center(cluster['point_list'], accuracy_fraction=80, min_acc=50)
        del cluster['point_list']
        cluster['time'] = '00:00:00'
        cluster['linecolor'] = '#ffffff'
        # cluster['trip'] = 0
        # cluster['drop'] = 0
        cluster['info_to_show_on_map'] = "np: " + str(cluster['np']) + ", nloc: " + str(cluster['nloc'])
        cluster['uidp'] = uidp

    # color by popularity
    maxlognloc = max([math.log(1 + float(cluster['nloc'])) for cluster in clusters])
    maxnloc = max([int(cluster['nloc']) for cluster in clusters])
    for cluster in clusters:
        cluster['markercolor'] = ut.range_color(math.log(1 + float(cluster['nloc'])) / maxlognloc)
        if maxnloc == cluster['nloc']:
            cluster['loctype'] = 'home'

    # sort by popularity
    clusters = sorted(clusters, key=operator.itemgetter('nloc', 'np', 'latitude', 'longitude'), reverse=True)

    # re-index
    new_idx = {}
    for idx, cluster in enumerate(clusters):
        cluster['id'] = idx
        new_idx[cluster['id_old']] = cluster['id']

    # refresh cluster index in the location data
    for loc in locs:
        loc['c_orig'] = new_idx[loc['c_orig']]

    # save clusters to file
    ut.dictlist2csv(clusters, PATHOUT + "location_clusters.csv")

    return str(num_clusters)


def process_structures(segs, point_dict, logger=None, silent=False):
    """ """
    # structure
    # struct_id - sid
    # coalesce locations first ->
    # coalesce trips first ->

    # tol_short_loc_time = 15, tol_short_loc_dist = 750,
    # tol_short_trip_time = 10, tol_short_trip_dist = 300,

    # pt_central_bangalore = (12.975158, 77.590022)

    """ initialize struct with segs """
    # pprint(segs)
    assert segs[-1]['type'] == 'end' or len(segs) == 1 and segs[-1]['type'] == 'gap'
    segs_struct = []
    for seg in segs[:-1]:
        seg_new = {
            'sid': seg['id'],
            'segs': [seg['id']],
            'trip': seg['type'],
            'date': seg['date'],
            'uidp': seg['uidp'],
            't00': seg['t00'],
            't01': seg['t01'],
            't00_': seg['t00_'],
            't01_': seg['t01_'],
            't_start_bg': seg['t0bg'],
            't_start_bg_': seg['t0bg_']
        }
        segs_struct.append(seg_new)

    # to navigate easily
    # access segments based on id
    # pprint(segs)
    # segs_dict = {seg['id']: seg for seg in segs if seg['type'] != 'end'}
    segs_dict = {seg['id']: seg for seg in segs}
    segs_next_id = {segs[idx]['id']: segs[idx+1]['id'] for idx in range(len(segs) - 1)}

    """ checks """
    # assert segs_struct[0]['type'] in ['loc', 'gap']

    """ START MERGING """
    idx_seg = 0
    while idx_seg < len(segs_struct) - 1:
        current_seg_struct = segs_struct[idx_seg]
        current_type = current_seg_struct['trip']

        next_seg_struct = segs_struct[idx_seg + 1]
        next_type = next_seg_struct['trip']
        next_dist = segs_dict[next_seg_struct['segs'][0]]['line_dist']
        next_dur = segs_dict[next_seg_struct['segs'][0]]['duration']

        if idx_seg + 2 < len(segs_struct):
            next2_seg_struct = segs_struct[idx_seg + 2]
            next2_type = next2_seg_struct['trip']
        else:
            next2_type = 'none'

        if current_type in ['loc', 'gap']:
            current_seg_struct['trip'] = 'loc'
            """ loc/gap + loc/gap """
            """ OR
                loc/gap + trip/jump/jumptrip + loc/gap if the middle trip is short
                or
                (next_type in ['trip', 'jump', 'jumptrip'] and
                next_dist < tol_short_trip_dist and next_dur < tol_short_trip_time and
                next2_type in ['none', 'gap', 'loc'])
            """
            if next_type in ['loc', 'gap']:
                next_seg = segs_dict[next_seg_struct['segs'][0]]
                if not silent:
                    logger.warning("Dropping LOC  | id=" + str(next_seg['id']).zfill(3) +
                                   " | " + next_seg['date'] + " " + next_seg['t01'] + "-" + next_seg['t11'] +
                                   " | start=" + str(next_seg['start']))
                # merge next into current
                # todo: eliminate the short trip conditon, it's too general
                current_seg_struct['segs'] += next_seg_struct['segs']
                del segs_struct[idx_seg + 1]
            else:
                idx_seg += 1
        else:
            assert current_type in ['trip', 'jump', 'jumptrip']
            current_seg_struct['trip'] = 'trip'
            """ trip + trip """
            """ OR """
            """ trip + loc/gap if loc if short"""
            """
                    or \
                    (next_type in ['loc', 'gap'] and
                     next_dur < tol_short_loc_time)
            """
            if next_type in ['trip', 'jump', 'jumptrip']:
                next_seg = segs_dict[next_seg_struct['segs'][0]]
                if not silent:
                    logger.warning("Dropping " + current_type.upper().replace("JUMPTRIP", "JUTR") +
                                   " | id=" + str(next_seg['id']).zfill(3) +
                                   " | " + next_seg['date'] + " " + next_seg['t01'] + "-" + next_seg['t11'] +
                                   " | start=" + str(next_seg['start']))
                # todo: remove location merging
                # merge next into current
                current_seg_struct['segs'] += next_seg_struct['segs']
                del segs_struct[idx_seg + 1]
            else:
                idx_seg += 1

    """ compute stats """
    for idx_seg, seg_struct in enumerate(segs_struct):
        seg_struct['types'] = ', '.join([segs_dict[idx]['type'] for idx in seg_struct['segs']])

        # pointlist
        seg_a = segs_dict[seg_struct['segs'][0]]
        idx_seg_z = seg_struct['segs'][-1]
        if idx_seg_z in segs_next_id:
            next_id = segs_next_id[idx_seg_z]
        else:
            max_id = max(segs_next_id.keys())
            assert idx_seg_z > max_id
            next_id = max_id

        """ location *diameter """
        # origin is the first point ON the segment.
        i_orig = seg_a['start']
        # destination is the last point on the segment (before the first point of the next segments)
        i_dest = min(segs_dict[next_id]['start'], max(point_dict))

        # for trips, include pt before, and pt after
        if seg_struct['trip'] == 'trip':
            if seg_a['type'] == 'trip':
                i_orig = get_prev_idx(i_orig, point_dict)
        else:
            # for locations, include only last point in the location
            i_dest = get_prev_idx(i_dest, point_dict)

        # save
        seg_struct['pts_od'] = [i_orig, i_dest]

        # diameter
        point_list_this_seg = [point_dict[idx] for idx in sorted(point_dict.keys()) if i_orig <= idx <= i_dest]
        points_latlon = [(float(point['latitude']), float(point['longitude'])) for point in point_list_this_seg]
        if len(points_latlon) > 1:
            seg_struct['center_diameter'], _ = diameter(points_latlon.copy())
        else:
            seg_struct['center_diameter'] = -1

        # for trips: straight line distance between the pt before and the pt after
        if seg_struct['trip'] == 'trip':
            seg_struct['line_dist'] = great_circle(points_latlon[0], points_latlon[-1]).meters

        """ precision """
        seg_struct['type_first'] = segs_dict[seg_struct['segs'][0]]['type']
        seg_struct['type_last'] = segs_dict[seg_struct['segs'][-1]]['type']

        """ TIMES """
        compute_times_struct(idx_seg=idx_seg, seg_struct=seg_struct,
                             segs_dict=segs_dict, segs_struct=segs_struct, point_dict=point_dict)

        """ for trips """
        if seg_struct['trip'] == 'trip':
            seg_struct['trip'] = 1
            # excludes first point of next segment
            seg_struct['pl'] = sum([segs_dict[idx]['path_length'] for idx in seg_struct['segs']
                                             if segs_dict[idx]['type'] in ['trip', 'jump', 'jumptrip']])

            seg_struct['plb'] = sum([segs_dict[idx]['path_length_inbang'] for idx in seg_struct['segs']
                                                    if segs_dict[idx]['type'] in ['trip', 'jump', 'jumptrip']])

            seg_struct['plb_jump'] = sum([segs_dict[idx]['path_length_inbang'] for idx in seg_struct['segs']
                                                    if segs_dict[idx]['type'] in ['jump', 'jumptrip']])

            # number of locations between 5 and 15 minutes
            seg_struct['n_loc'] = 0
            for idx in seg_struct['segs']:
                if segs_dict[idx]['type'] in ['loc', 'gap'] and 5 <= segs_dict[idx]['duration']:
                    seg_struct['n_loc'] += 1

        # """ for locations """
        else:
            assert seg_struct['trip'] == 'loc'
            seg_struct['trip'] = 0
            # get location structure center
            seg_struct['center_lat'], seg_struct['center_lon'], seg_struct['center_acc'], seg_struct['center_std'] = \
                pc.point_list_center(point_list_this_seg, accuracy_fraction=80, min_acc=0)

            pt_latlong = (seg_struct['center_lat'], seg_struct['center_lon'])
            seg_struct['center_d2c'] = great_circle(pt_latlong, pt_central_bangalore).meters

    # trip location
    for idx, seg_struct in enumerate(segs_struct):
        if seg_struct['trip'] == 1 and 0 < idx:
            assert segs_struct[idx - 1]['trip'] == 0
            seg_struct['center_orig_lat'], seg_struct['center_orig_lon'] = segs_struct[idx - 1]['center_lat'], \
                                                                           segs_struct[idx - 1]['center_lon']

            # distance between origin/destination and central Bangalore
            seg_struct['center_orig_d2b'] = great_circle((seg_struct['center_orig_lat'], seg_struct['center_orig_lon']),
                                                         (pt_central_bangalore)).meters

        if seg_struct['trip'] == 1 and idx < len(segs_struct) - 1:
            assert segs_struct[idx + 1]['trip'] == 0
            seg_struct['center_dest_lat'], seg_struct['center_dest_lon'] = segs_struct[idx + 1]['center_lat'], \
                                                                           segs_struct[idx + 1]['center_lon']

            # distance between origin/destination and central Bangalore
            seg_struct['center_dest_d2b'] = great_circle((seg_struct['center_dest_lat'], seg_struct['center_dest_lon']),
                                                         (pt_central_bangalore)).meters

        if seg_struct['trip'] == 1 and 0 < idx < len(segs_struct) - 1:
            # direct straight line distance between orig center and dest center
            seg_struct['center_line_dist'] = great_circle(
                (seg_struct['center_orig_lat'], seg_struct['center_orig_lon']),
                (seg_struct['center_dest_lat'], seg_struct['center_dest_lon'])).meters

    return segs_struct


def find_chain_end(idx_first, segs_struct, locations, max_loc_dist, max_loc_dur_chain):

    # first segment
    seg = segs_struct[idx_first]
    assert seg['trip'] == 1
    if idx_first == 0:
        chain_orig = (None, None)
    else:
        chain_orig = (seg['center_orig_lat'], seg['center_orig_lon'])

    # common locations
    home_pt, work_pt, work2_pt, work3_pt = locations

    # start going through segments
    idx_current = idx_first
    chain_stop = False
    chain_properties = {}
    while not chain_stop:
        # if one of the stop conditions satisfied, stop here:
        # - NEXT location is > 30 minutes
        # - NEXT location is known location
        # - NEXT location is same as origin

        # decide to stop during the last trip
        if segs_struct[idx_current]['trip'] == 1:
            loc_latlong = (segs_struct[idx_current]['center_dest_lat'], segs_struct[idx_current]['center_dest_lon'])

            # close to a make_hw_maps location?
            dist2h = great_circle(home_pt, loc_latlong).meters < max_loc_dist
            dist2w = great_circle(work_pt, loc_latlong).meters < max_loc_dist
            dist2w2 = dist2w3 = False
            if work2_pt: dist2w2 = great_circle(work2_pt, loc_latlong).meters < max_loc_dist
            if work3_pt: dist2w3 = great_circle(work3_pt, loc_latlong).meters < max_loc_dist

            # is last or next is last?
            c0_last = idx_current == (len(segs_struct) - 1)
            if idx_current == (len(segs_struct) - 2) and segs_struct[-1]['trip'] == 0: c0_last = True

            # is the next location long? (or is it the last?)
            if c0_last: c1_next_long_loc = False
            else: c1_next_long_loc = segs_struct[idx_current + 1]['dur_bg_mm'] > max_loc_dur_chain
            # is the trip chain destination a make_hw_maps location?
            c2_next_close_main_loc = dist2h or dist2w or dist2w2 or dist2w3
            # is the trip chain origin same as the destination?
            c3_next_isorigin = chain_orig[0] and (great_circle(chain_orig, loc_latlong).meters < max_loc_dist)
            # print("idx " + str(idx_current) + " dest " + str(loc_latlong) +
            #       " origin " + str(chain_orig) + " dist " + str(great_circle(chain_orig, loc_latlong).meters))

            if c0_last or c1_next_long_loc or c2_next_close_main_loc or c3_next_isorigin:
                chain_stop = True
                chain_properties = {'dest_main_loc': c2_next_close_main_loc, 'orig_is_dest': c3_next_isorigin}
            else:
                idx_current += 1
        else:
            # skip locations, except if last:
            # c0_last = idx_current == (len(segs_struct) - 1)
            # if c0_last:
            #     chain_stop = True
            #     chain_properties = {'dest_main_loc': c2_next_close_main_loc, 'orig_is_dest': c3_next_isorigin}
            # else:
            idx_current += 1
    return idx_current, chain_properties


def process_chain(segs_struct, locations, point_dict, max_loc_dur_mm=30, max_loc_dist=300, logger=None, silent=True):
    """
    max_loc_dur_chain = 30 minutes max
    max_loc_dist = 300 meters maximum distance to make_hw_maps location
    """

    """ Go through and tag chained trips """
    idx_chain = 0
    idx_struct = 0
    segs_chain = []

    while idx_struct < len(segs_struct):
        if not silent:
            logger.info("processing idx_struct " + str(idx_struct))
        seg = segs_struct[idx_struct]

        assert seg['trip'] in [0, 1]

        """ Location -> add as non-chain """
        if seg['trip'] == 0:
            seg_new = dict()
            for field in ['trip', 'date', 'uidp', 'dur_bg_mm', 't_start_0', 't_start_1', 't_stop_0', 't_stop_1',
                          't_start_bg', 't_start_bg_', 't_stop_bg', 't_stop_bg_']:
                seg_new[field] = seg[field]
            seg_new['sids'] = [seg['sid']]
            seg_new['segs'] = seg['segs']
            seg_new['chain'] = idx_chain
            seg_new['chain_diameter'] = int(np.round(float(seg['center_diameter'])))
            segs_chain.append(seg_new)

            # increment
            idx_struct += 1
            # increment chain
            idx_chain += 1

        """ Trip -> start and follow and complete new chain """
        if seg['trip'] == 1:
            idx_last, chain_properties = find_chain_end(idx_first=idx_struct, segs_struct=segs_struct,
                                                        locations=locations,
                                                        max_loc_dist=max_loc_dist,
                                                        max_loc_dur_chain=max_loc_dur_mm)
            #
            seg_new = dict()
            for field in ['trip', 'date', 'uidp', 't_start_0', 't_start_1', 't_start_bg', 't_start_bg_']:
                seg_new[field] = seg[field]
            for field in ['t_stop_0', 't_stop_1', 't_stop_bg', 't_stop_bg_']:
                seg_new[field] = segs_struct[idx_last][field]

            # total length and duration (for duration: including locations)
            for field in ['pl', 'plb', 'dur_bg_mm']:
                seg_new[field] = sum([segs_struct[idx].get(field, 0) for idx in range(idx_struct, idx_last + 1)])

            # define total duration *on trips*
            seg_new['dur_trips'] = sum([segs_struct[idx].get('dur_bg_mm', 0)
                                        for idx in range(idx_struct, idx_last + 1)
                                        if segs_struct[idx]['trip'] == 1])

            try:
                seg_new['center_orig_d2b'] = seg.get('center_orig_d2b', 0)
                seg_new['center_orig_lat'] = seg.get('center_orig_lat', '')
                seg_new['center_orig_lon'] = seg.get('center_orig_lon', '')
                seg_new['center_dest_d2b'] = segs_struct[idx_last].get('center_dest_d2b', 0)
                seg_new['center_dest_lat'] = segs_struct[idx_last].get('center_dest_lat', '')
                seg_new['center_dest_lon'] = segs_struct[idx_last].get('center_dest_lon', '')

                # linear distance between endpoints
                if seg_new['center_orig_lat'] and seg_new['center_dest_lat']:
                    seg_new['center_line_dist'] = great_circle(
                        (seg_new['center_orig_lat'], seg_new['center_orig_lon']),
                        (seg_new['center_dest_lat'], seg_new['center_dest_lon'])).meters

                # entire diameter
                if idx_struct == idx_last:
                    seg_new['chain_diameter'] = int(np.round(float(seg['center_diameter'])))
                else:
                    [i_orig, _] = eval(str(seg['pts_od']))
                    [_, i_dest] = eval(str(segs_struct[idx_last]['pts_od']))

                    # diameter
                    point_list_this_seg = [point_dict[idx] for idx in sorted(point_dict.keys()) if
                                           i_orig <= idx <= i_dest]
                    points_latlon = [(float(point['latitude']), float(point['longitude'])) for point in
                                     point_list_this_seg]
                    if len(points_latlon) > 1:
                        seg_new['chain_diameter'], _ = diameter(points_latlon.copy())
                        seg_new['chain_diameter'] = int(np.round(seg_new['chain_diameter']))
                    else:
                        seg_new['chain_diameter'] = -1

            except Exception as e:
                print(seg['uidp'])
                raise e

            seg_new['sids'] = [segs_struct[idx]['sid'] for idx in range(idx_struct, idx_last + 1)]
            seg_new['segs'] = [seg_id for idx in range(idx_struct, idx_last + 1) for seg_id in segs_struct[idx]['segs']]
            seg_new['n_struct_loc'] = len([1 for idx in range(idx_struct, idx_last + 1)
                                           if segs_struct[idx]['trip'] == 0])
            seg_new['chain'] = idx_chain
            idx_chain += 1
            seg_new = {**seg_new, **chain_properties}
            segs_chain.append(seg_new)
            if not silent:
                logger.warning("starting chain " + str(idx_chain) + " struct " + str(idx_struct))

            # increment
            idx_struct = idx_last + 1

    return segs_chain


def process_segments_p1(deviceid, uidp, name, points_by_date, prefix="", logger=None, silent=False):
    # read dates
    sorted_dates = sorted([date for date in list(points_by_date.keys())])

    # keep track of unreasonable speeds in this dictionary
    distances_removed = []

    allsegs = {}
    points_list_active = {}
    points_dict = {}
    points_dict_active = {}

    for mydate in sorted_dates:
        prefix_date = prefix + mydate + " "
        if not silent: logger.info(prefix_date)
        allpoints = points_by_date[mydate]

        # continue only for non-dropped points
        assert set([point.get('drop', '0') for point in allpoints]).issubset({'', '0', '1'})
        points_list_active[mydate] = [point for point in allpoints if point.get('drop', '0') != '1']

        # dictionaries
        points_dict[mydate] = {int(line['idx']): line for line in allpoints}
        points_dict_active[mydate] = {int(line['idx']): line for line in points_list_active[mydate]}

        if not points_list_active[mydate]:
            # none
            segs = []

        else:
            # start with location always
            points_list_active[mydate][0]["trip"] = "0"

            # add token point at day endpoints
            point_0000 = points_list_active[mydate][0].copy()
            point_0000["time"] = "00:00:00"
            point_0000["trip"] = "0"
            point_0000["idx"] = 0

            point_1159 = points_list_active[mydate][-1].copy()
            point_1159['time'] = "23:59:59"
            point_1159["idx"] = int(points_list_active[mydate][-1]['idx']) + 1

            # add to list and dictionary
            points_list_active[mydate] = [point_0000] + points_list_active[mydate] + [point_1159]
            points_dict_active[mydate][0] = point_0000
            points_dict_active[mydate][point_1159["idx"]] = point_1159

            """ refresh dist/time 2prev/2next """
            point_simple_stats(points_list_active[mydate], points_dict_active[mydate])
            if not silent: logger.info(prefix_date + 'SEGMENT CODE: simple stats - done')

            """ segmentize - naive """
            segs = find_segs_naive(points_list_active[mydate])
            if not silent: logger.info(prefix_date + 'SEGMENT CODE: naive segments - done')

            """ add previously dropped points within location in points_list_active_dict[mydate] """
            process_ok_low_acc(segs, points_dict_active[mydate], points_dict[mydate], logger=logger, silent=silent)
            if not silent: logger.info(prefix_date + 'SEGMENT CODE: adding low acc in range of loc - done')

            # update the point list
            points_list_active[mydate] = [points_dict_active[mydate][idx]
                                          for idx in sorted(points_dict_active[mydate].keys())]

            """ update dist2next etc including low acc points """
            point_simple_stats(points_list_active[mydate], points_dict_active[mydate], distances_removed)

            """ detect gaps and jumps """
            segs = process_segs_gaps(segs, points_dict_active[mydate], logger=logger, silent=silent)
            if not silent: logger.info(prefix_date + 'SEGMENT CODE: detect gaps and jumps - done')

            # calculate location stats
            segs_stats(segs, points_dict_active[mydate])
            if not silent: logger.info(prefix_date + 'SEGMENT CODE: segment stats - done')

            # eliminate short trips
            if segs:
                if not silent: logger.info(prefix_date + 'SEGMENT CODE: segment drop short trips - len = ' + str(len(segs)))
                segs = process_short_trips(segs, points_dict_active[mydate], logger=logger, silent=silent)
                if not silent: logger.info(prefix_date + 'SEGMENT CODE: segment drop short trips - done - len = ' + str(len(segs)))

            """ ADDING IDS - done only once """
            segs_stats_ids(segs)

            # id
            for seg in segs:
                seg['name'] = name
                seg['deviceid'] = deviceid
                seg['uidp'] = uidp

        allsegs[mydate] = segs

    return allsegs, points_dict_active, sorted_dates, distances_removed


def process_segments_p2(ROOT_PATH, uidp, allsegs, points_dict, sorted_dates, chains_max_loc_dur_mm=30,
                        PATHIN='', PATHOUT='', logger=None, silent=False):
    """ """
    if not silent: logger.debug("Part 2 : processing : " + uidp)
    """ Read home and work(s) destinations from file """

    # read confirmed coded hw locations
    confirmed_hw_file = ROOT_PATH + 'data/raw_gps_homework/hw_confirmed.csv'
    confirmed_hw = ut.csv2dict(confirmed_hw_file, key_name='uidp')
    confirmed_hw = confirmed_hw.get(uidp, None)

    # read common locations
    locations_file = PATHIN + 'location_clusters.csv'
    locations = ut.csv2dict(locations_file)
    far_away_pt = {'latitude': 0, 'longitude': 0}
    if len(locations) > 0: loc0 = locations[0]
    else: loc0 = far_away_pt
    if len(locations) > 1: loc1 = locations[1]
    else: loc1 = far_away_pt

    # H and W locations
    home_pt = float(loc0['latitude']), float(loc0['longitude'])
    work_pt = float(loc1['latitude']), float(loc1['longitude'])
    work2_pt = None
    work3_pt = None

    # override common locations
    if confirmed_hw:
        if confirmed_hw['orig']:
            home_pt = eval(confirmed_hw['orig'])
        if confirmed_hw['dest']:
            work_pt = eval(confirmed_hw['dest'])
        if confirmed_hw['dest2']:
            work2_pt = eval(confirmed_hw['dest2'])
        if confirmed_hw['dest3']:
            work3_pt = eval(confirmed_hw['dest3'])
    else:
        if not silent: print("using common locations HW")
    locations = home_pt, work_pt, work2_pt, work3_pt

    """ Struct """
    allsegs_struct = []
    allsegs_chain = []
    for mydate in sorted_dates:
        segs = allsegs[mydate]

        # combine segments into trip and location structures
        if segs:
            segs_struct = process_structures(segs, points_dict[mydate], logger=logger, silent=silent)
            allsegs_struct += segs_struct

            # process chained trips
            segs_chain = process_chain(segs_struct, max_loc_dur_mm=chains_max_loc_dur_mm,
                                       point_dict=points_dict[mydate],
                                       locations=locations, logger=logger, silent=silent)
            allsegs_chain += segs_chain

    """ Print structures in a single file """
    if ut.dictlist2csv(allsegs_struct, PATHOUT + 'segs_struct.csv') and not silent:
        logger.info("written to file successfully " + PATHOUT + 'segs_struct.csv')

    # print all structures (PRETTY)
    if pretty_print_struct(allsegs_struct, PATHOUT + 'segs_struct.txt') and not silent:
        logger.info("written to file successfully " + PATHOUT + 'segs_struct.txt (PRETTY)')

    # print all chains
    if ut.dictlist2csv(allsegs_chain, PATHOUT + 'segs_chain_' + str(chains_max_loc_dur_mm) + '.csv') and not silent:
        logger.info("written to file successfully " + PATHOUT + 'segs_chain.csv')

    # print all segments
    allsegs = [line for mydate in sorted(allsegs.keys()) for line in allsegs[mydate]]
    if ut.dictlist2csv(allsegs, PATHOUT + 'trips.csv') and not silent:
        logger.info("written to file successfully " + PATHOUT + 'trips.csv')

    return True


def tp_segments(deviceid, uidp, name,
                chains_max_loc_dur_mm=30,
                ROOT_PATH=None, PATHIN_STUB=None, PATHOUT_STUB=None,
                prefix="", dates=None, debug=False, silent=False):
    """
    By date, read points and:
    --- classify naively into location/trip segments
    --- identify gaps, jumps, tripjumps
    --- eliminate short trips

    --- find clusters of locations
    """
    if silent:
        logger = None
        print(prefix + "processing uidp " + uidp)
    else:
        logger = logging.getLogger()
        now = str(datetime.datetime.now())[0:19].replace(" ", "_").replace(":", "")
        log_file_now = ROOT_PATH + 'data/logs/seg_struct/tp_segments_struct_' + now + '.txt'
        ul.enable_logging(logger, log_file=log_file_now)
        logger.info("Starting log file for TP SEGMENTS. DEBUG = " + str(debug))

    if not prefix:
        prefix = uidp + " "

    if not silent: logger.info('Segment coding for deviceid = ' + deviceid + ', uidp = ' + uidp + ', ' + name)

    """ paths """
    PATHIN_FOLDER = PATHIN_STUB + uidp + "/"
    PATHOUT_FOLDER = PATHOUT_STUB + uidp + "/"
    assert PATHOUT_FOLDER
    if debug:
        PATHOUT_FOLDER = ROOT_PATH + "data/sandbox_trips/"
    if not os.path.exists(PATHOUT_FOLDER + 'trips/'):
        os.makedirs(PATHOUT_FOLDER + 'trips/')

    # read all points for this period
    if not silent: logger.info(prefix + PATHIN_FOLDER + 'points/')
    points_by_date = ut.csv2pointlist(PATHIN_FOLDER + 'points/', silent=silent)

    if not points_by_date:
        if not silent: logger.info("No points data. Done. " + uidp)
        return []
    else:
        if dates:
            points_by_date = {mydate: points_by_date[mydate] for mydate in points_by_date
                              if mydate in dates}
            if not silent: logger.info("Restricted to dates in dates")

        if not silent:
            logger.info(prefix + "SEGMENT CODE: read points - done")
            logger.info('Number of dates with point data = ' + str(len(points_by_date)))

        """ Part 1: naive, gaps, eliminate short trips """
        allsegs, points_dict_active, sorted_dates, distances_removed = \
            process_segments_p1(deviceid, uidp, name, points_by_date, prefix=prefix, logger=logger, silent=silent)
        if not silent: logger.info('Naive coding yields ' + str(len(allsegs)) + 'segments.')

        """ Part 2: detect clusters """
        alllocs = [seg for mydate in allsegs for seg in allsegs[mydate] if seg['type'] == 'loc']
        if not silent: logger.info(prefix + " >>> detecting clusters <<<")
        if alllocs:
            allsegs_list = [seg for mydate in allsegs for seg in allsegs[mydate]]
            # detect make_hw_maps locations - cluster
            n_clusters = process_loc_clusters(uidp, allsegs_list, PATHOUT=PATHIN_FOLDER, logger=logger, silent=silent)
            if not silent: logger.info('Clusters ' + str(n_clusters))

            segs_stats_cluster_codes(allsegs_list)

        """ Part 3: meta trips, tag OK trips, c """
        if not silent: logger.info(prefix + " >>> segments part 2 <<<")
        process_segments_p2(ROOT_PATH, uidp, allsegs, points_dict_active, sorted_dates,
                            chains_max_loc_dur_mm=chains_max_loc_dur_mm,
                            PATHIN=PATHIN_FOLDER, PATHOUT=PATHOUT_FOLDER, logger=logger, silent=silent)

        now_str = str(datetime.datetime.now())[0:19]
        if silent:
            print("uidp " + uidp + " finished at " + now_str)
        else:
            logger.info("uidp " + uidp + " finished at " + now_str)
        return distances_removed


def tp_chains(ROOT_PATH, uidp, max_loc_dur_mm, PATHIN, PATHOUT=None, struct_file=None, debug=False, silent=True):
    """
    max_loc_dur_mm = maximum duration at a certain location to be included in a chain
    (in minutes)
    """

    # paths and files
    if not PATHOUT:
        PATHOUT = PATHIN + uidp + '/'
    if not struct_file:
        struct_file = PATHIN + uidp + '/segs_struct.csv'

    """ Logger """
    logger = logging.getLogger()
    now = str(datetime.datetime.now())[0:19].replace(" ", "_").replace(":", "")
    log_file_now = ROOT_PATH + 'data/logs/seg_struct/tp_segments_chains_' + now + '.txt'
    ul.enable_logging(logger, log_file=log_file_now)
    if not silent:
        logger.info("Starting log file for TP CHAINS. DEBUG = " + str(debug))

    """ Read points into dict dict """
    points_by_date = ut.csv2pointlist(PATHIN + uidp + '/points/', silent=True)
    for mydate in points_by_date:
        points_by_date[mydate] = {int(line['idx']): line for line in points_by_date[mydate]}

    """ Read struct trips """
    segs_struct = ut.csv2dict(struct_file)

    # list of segments
    for seg in segs_struct:
        seg['segs'] = eval(seg['segs'])
        seg['trip'] = int(seg['trip'])
        seg['sid'] = int(seg['sid'])
        for field in ['pl', 'plb', 'dur_bg_mm']:
            if field in seg:
                if seg[field] == '': del seg[field]  # seg['field'] = 0
                else: seg[field] = float(seg[field])

        for field in ['center_orig_lat', 'center_orig_lon', 'center_dest_lat', 'center_dest_lon']:
            if field in seg:
                if seg[field] == '': del seg[field]
                else: seg[field] = float(seg[field])

    """ Read home and work(s) destinations from file """
    # read confirmed coded hw locations
    confirmed_hw_file = ROOT_PATH + 'data/raw_gps_homework/hw_confirmed.csv'
    confirmed_hw = ut.csv2dict(confirmed_hw_file, key_name='uidp')
    confirmed_hw = confirmed_hw.get(uidp, None)

    # read common locations
    locations_file = PATHIN + 'location_clusters.csv'
    locations = ut.csv2dict(locations_file)
    far_away_pt = {'latitude': 0, 'longitude': 0}
    if len(locations) > 0:
        loc0 = locations[0]
    else:
        loc0 = far_away_pt
    if len(locations) > 1:
        loc1 = locations[1]
    else:
        loc1 = far_away_pt

    # H and W locations
    home_pt = float(loc0['latitude']), float(loc0['longitude'])
    work_pt = float(loc1['latitude']), float(loc1['longitude'])
    work2_pt = None
    work3_pt = None

    # override common locations
    if confirmed_hw:
        if confirmed_hw['orig']:
            home_pt = eval(confirmed_hw['orig'])
        if confirmed_hw['dest']:
            work_pt = eval(confirmed_hw['dest'])
        if confirmed_hw['dest2']:
            work2_pt = eval(confirmed_hw['dest2'])
        if confirmed_hw['dest3']:
            work3_pt = eval(confirmed_hw['dest3'])
    else:
        print("using common locations HW")
    locations = home_pt, work_pt, work2_pt, work3_pt

    """ CHAINS """
    allsegs_chain = []
    sorted_dates = sorted(set([seg['date'] for seg in segs_struct]))
    for mydate in sorted_dates:
        # print("Processing " + mydate)
        segs = [seg for seg in segs_struct if seg['date'] == mydate]

        # combine segments into trip and location structures
        if segs:
            # process chained trips
            segs_chain = process_chain(segs, locations=locations,
                                       point_dict=points_by_date[mydate],
                                       max_loc_dur_mm=max_loc_dur_mm, max_loc_dist=300,
                                       logger=logger)
            allsegs_chain += segs_chain

    # print all chains
    if ut.dictlist2csv(allsegs_chain, PATHOUT + 'segs_chain_' + str(max_loc_dur_mm) + '.csv'):
        if not silent:
            logger.info("written to file successfully " + PATHOUT + 'segs_chain.csv')

    return True


def collect_segs(ROOT_PATH, PATH_OUT, PATHIN_STUB=None, uidps=None, refresh_meta=True, collect_chains=None):
    """
    Simple script to read individual trip files (by uidp and date)
    collect_chains should contain the number (15, 30 or 60) that identifies the
    """

    PATH_UIDP_DEVICEID = ROOT_PATH + "data/coded_cto/crosswalk_uidp_deviceid.csv"
    crosswalk_dict = ut.csv2dict(PATH_UIDP_DEVICEID)

    if uidps:
        crosswalk_dict = [subject for subject in crosswalk_dict if subject['uidp'] in uidps]

    clusters = []
    allsegs_meta = []
    allsegs_struct = []
    allsegs_simple = []
    allsegs_chains = []

    if not PATHIN_STUB:
        PATHIN_STUB = ROOT_PATH + "data/coded_gps_byid/"

    for subject in crosswalk_dict:
        deviceid, uidp, name = subject["deviceid"], subject["uidp"], subject["name"]
        print("COLLECT: Processing segments " + str([deviceid, uidp, name]))

        ID_PATH = PATHIN_STUB + uidp + "/"

        # read clusters
        clusters += ut.csv2dict(ID_PATH + 'location_clusters.csv')

        # read meta trips (glob them or just read the overall)
        meta_file = ID_PATH + 'trips_meta.csv'
        allsegs_meta_this = ut.csv2dict(meta_file)

        if refresh_meta:
            PATH_READ_META = glob.glob(ID_PATH + "trips/trip_meta_" + uidp + "_*.csv")
            allsegs_meta_this = []
            for PATH_FILE in PATH_READ_META:
                segs = ut.csv2dict(PATH_FILE)
                # add unique ID within segs
                for idx, seg in enumerate(segs):
                    seg['id_new'] = idx
                allsegs_meta += segs
                allsegs_meta_this += segs

            # TRIPS (SIMPLE)
        simple_file = ID_PATH + 'trips.csv'
        allsegs_simple += ut.csv2dict(simple_file)

        # TRIPS (STRUCT)
        struct_file = ID_PATH + 'segs_struct.csv'
        allsegs_struct += ut.csv2dict(struct_file)

        # CHAINS
        if collect_chains:
            chains_file = ID_PATH + 'segs_chain_' + str(collect_chains) + '.csv'
            allsegs_chains += ut.csv2dict(chains_file)

        # write THIS in the UIDP folder
        ut.dictlist2csv(allsegs_meta_this, ID_PATH + 'trips_meta.csv')

    # write ALL in the global folder
    ut.dictlist2csv(clusters, PATH_OUT + 'clusters.csv')
    ut.dictlist2csv(allsegs_meta, PATH_OUT + 'trips_meta.csv')
    ut.dictlist2csv(allsegs_struct, PATH_OUT + 'trips_struct.csv')
    ut.dictlist2csv(allsegs_simple, PATH_OUT + 'trips.csv')
    ut.dictlist2csv(allsegs_chains, PATH_OUT + 'chains_' + str(collect_chains) + '.csv')
    return True



def time_overlap(myseg, t0, t1):
    """
    Compute the overlap between the segments t0,t1 and s0,s1 and return that amount if it's trip, loc or gap<60minutes
    :param myseg:
    :param t0: - fractional hour
    :param t1: - fractional hour
    :return:
    """
    s0 = ut.hf_hhmm(myseg['t01'])
    s1 = ut.hf_hhmm(myseg['t11'])
    # s0 = ut.hf_hhmmss(myseg['t_start_bg'])
    # s1 = ut.hf_hhmmss(myseg['t_stop_bg'])

    assert s0 <= s1
    assert t0 <= t1

    if s1 <= t0:
        olap = 0.0
    elif s0 <= t0 <= s1:
        rend = min(s1, t1)
        olap = rend - t0
    elif s0 <= s1 <= t1:
        olap = s1 - s0
    elif s0 <= t1 < s1:
        olap = t1 - s0
    else:
        assert t1 < s0
        olap = 0.0

    assert olap <= t1 - t0

    return olap


def flag_outstation(segs, debug=False):
    """
    Outstation definition: 50% of time outstation
    :param segs:
    :return:
    """
    # start and end day outside Bangalore
    start_outstation = False
    end_outstation = False
    locs = [seg for seg in segs if seg['type'] == 'loc']
    # trips = [seg for seg in segs if seg['type'] == 'trip']
    # trips_ok = [seg for seg in segs if seg['trip_ok'] == '1']
    if locs:
        start_outstation = int(locs[0]['outstation_full']) == 1
        end_outstation = int(locs[-1]['outstation_full']) == 1

    time_out = sum([time_overlap(seg, 7, 21) for seg in segs if int(seg['outstation_full']) == 1])
    time_tot = sum([time_overlap(seg, 7, 21) for seg in segs])
    if time_tot == 0.0:
        return False
    elif time_out / time_tot > 0.5 and (start_outstation or end_outstation):
        print("returning true")
        return True
    else:
        return False


def gap_type(gap):
    dur, dist = float(gap['duration']), float(gap['line_dist'])

    if 1.5 * 60 <= dur < 4 * 60 and   0 <= dist < 400:
        return 'med'
    if 0.5 * 60 <= dur < 2 * 60 and 400 <= dist < 2000:
        return 'med'
    if 2.0 * 60 <= dur          and 400 <= dist < 2000:
        return 'gap'
    if 4.0 * 60 <= dur          and   0 <= dist < 2000:
        return 'gap'
    if 2000 <= dist:
        return 'jump'

    return 'low'


def compute_qm(uidp, segs_today, mydate):
    """
    Quality measures for data for a single person
    --- based on the meta trip data
    --- segment penalty for each segment (adds up to at most 2.0 or 2.5 for absolutely no data
    --- also other measures (older)
    """

    segs_inrange = [seg for seg in segs_today if time_overlap(seg, 7, 21) > 0.0]

    qm_dict = {
        'uidp': uidp,
        'date': mydate,
        'gap_hr': 0.0 - (not segs_inrange),
        'gap_n': 0.0 - (not segs_inrange),
        'gap_dist': 0.0 - (not segs_inrange),
        'jump_hr': 0.0 - (not segs_inrange),
        'jump_dist': 0.0 - (not segs_inrange),
        'outstation': int(flag_outstation(segs_today))
    }

    # determine if they left from home and returned home
    for idx, myseg in enumerate(segs_today):
        olap = time_overlap(myseg, 7, 21)

        if myseg['type'] == 'gap':
            olap = max(0.0, olap - 0.75)
            qm_dict['gap_hr'] += olap

        # jumptrip with degenerate trip is "pre"
        if myseg['type'] == 'jumptrip':
            qm_dict['jump_hr'] += olap
            try:
                if float(segs_today[idx - 1]['line_dist']) < 100:
                    qm_dict['jump_dist'] += float(myseg['line_dist']) / 1000.0
            except Exception as e:
                raise e

        if myseg['type'] == 'jump':
            qm_dict['jump_hr'] += olap
            qm_dict['jump_dist'] += float(myseg['line_dist']) / 1000.0

    if qm_dict['gap_hr'] == -1.0:
        qm_dict['tot_gap_hr'] = -1.0
        qm_dict['data_quality'] = 'no data'
    else:
        # qm_dict['tot_gap_hr'] = qm_dict['gap_hr'] + qm_dict['jump_pre_hr'] + qm_dict['jump_post_hr']
        qm_dict['tot_gap_hr'] = qm_dict['gap_hr'] + qm_dict['jump_hr']
        # 2.5 hours gap OR 1.5 km jump pre
        # 3.5 hours *total* gap OR 1.5 KM jump pre
        if qm_dict['tot_gap_hr'] > 3.5 or qm_dict['jump_dist'] > 1.5:
            qm_dict['data_quality'] = 'bad'
            # 5.5 hours gap
            # 6.5 hours *total* gap
            if qm_dict['tot_gap_hr'] > 6.5:
                qm_dict['data_quality'] = 'halfday'
        else:
            qm_dict['data_quality'] = 'good'

    return qm_dict
