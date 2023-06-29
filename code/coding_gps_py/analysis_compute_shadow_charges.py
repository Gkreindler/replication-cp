import datetime
from operator import itemgetter
import utils as ut
import utils_geometry as utg
import time
from pprint import pprint
from math import floor

from geopy.distance import great_circle


def seg_dist(seg):
    """
    Return best estimate of segment total distance
    :param seg:
    :return:
    """

    return float(seg.get('path_length_inbang_tot', 0.0) or 0.0) / 1000

    # if seg.get('path_length_trips', 0.0):
    #     assert seg['path_length'] not in ['', '0', '0.0']
    #     return float(seg.get('path_length_trips', 0.0) or 0.0) / 1000
    #
    # return max(float(seg['line_dist'] or 0.0),
    #            float(seg.get('path_length', '0.0') or 0.0)) / 1000


def pprint_segs(segs):
    print(('type', 'start', 'end', 'dist'))
    for seg in segs:
        print((seg['type'], seg['t01'], seg['t11'], "{:3.1f}".format(seg_dist(seg))))


def time_overlap(myseg, t0, t1):
    """
    Compute the overlap between the segments t0,t1 and s0,s1 and return that amount if it's trip, loc or gap<60minutes
    :param myseg:
    :param t0:
    :param t1:
    :return:
    """
    s0 = ut.hf_hhmm(myseg['t01'])
    s1 = ut.hf_hhmm(myseg['t11']) + 1 / 60
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


def flag_no_data(segs, tam0, tpm0):
    # compute data overlap with AM charges (+/-1h and ramp only, resp)
    hr_missing = 0.0 + 10 * (not segs)
    hr_missing_am = 0.0 + 10 * (not segs)
    hr_missing_pm = 0.0 + 10 * (not segs)

    for myseg in segs:
        olap = time_overlap(myseg, 7, 21)
        olap_am = time_overlap(myseg, 7, 13)
        olap_pm = time_overlap(myseg, 15, 21)

        mytype = myseg['type']
        if mytype == 'gap':
            olap = max(0.0, olap - 0.75)
            olap_am = max(0.0, olap_am - 0.75)
            olap_pm = max(0.0, olap_pm - 0.75)

        if olap > 0.0 and mytype in ['gap', 'jumptrip', 'jump'] and myseg.get('gap_ok', '') != '1':
            hr_missing += olap
            hr_missing_am += olap_am
            hr_missing_pm += olap_pm

    # good_overlap_am = sum([good_time_overlap(seg, tam0 - 1, tam0 + 4) for seg in segs])
    # good_overlap_ramp_am = sum([good_time_overlap(seg, tam0, tam0 + 3) for seg in segs])
    total_jumps_am = sum([seg_dist(seg) for seg in segs
                          if time_overlap(seg, tam0, tam0 + 3) and seg['type'] in ['jump', 'jumptrip']])
    # # if 60 minute missing from 5h or 45 minutes missing form 3h
    # flag_missing_data_am = good_overlap_am < 5 - 1 or \
    #                        good_overlap_ramp_am < 3 - 0.75 or \
    #                        total_jumps_am >= 5000
    #
    # # compute data overlap with PM charges (+/-1h and ramp only, resp)
    # good_overlap_pm = sum([good_time_overlap(seg, tpm0 - 1, tpm0 + 4) for seg in segs])
    # good_overlap_ramp_pm = sum([good_time_overlap(seg, tpm0, tpm0 + 3) for seg in segs])
    total_jumps_pm = sum([seg_dist(seg) for seg in segs
                          if time_overlap(seg, tpm0, tpm0 + 3) and seg['type'] in ['jump', 'jumptrip']])
    #
    # # if 60 minute missing from 5h or 45 minutes missing form 3h
    # flag_missing_data_pm = good_overlap_pm < 5 - 1 or \
    #                        good_overlap_ramp_pm < 3 - 0.75 or \
    #                        total_jumps_pm >= 5000

    # print(hr_missing_am)
    # print(total_jumps_am)

    good_overlap_am = hr_missing_am < 2.0
    good_overlap_pm = hr_missing_pm < 2.0

    flag_missing_data_am = (not good_overlap_am) or (total_jumps_am > 4.0)
    flag_missing_data_pm = (not good_overlap_pm) or (total_jumps_pm > 4.0)

    stats = {
        'flag_missing_data_am': flag_missing_data_am,
        'flag_missing_data_pm': flag_missing_data_pm,
        'good_overlap_am': good_overlap_am,
        # 'good_overlap_ramp_am': good_overlap_ramp_am,
        'total_jumps_am': total_jumps_am,
        'good_overlap_pm': good_overlap_pm,
        # 'good_overlap_ramp_pm': good_overlap_ramp_pm,
        'total_jumps_pm': total_jumps_pm,
    }
    # pprint(stats)

    return flag_missing_data_am, flag_missing_data_pm, stats


def flag_no_trip(segs, tam_0, tam_1, tpm_0, tpm_1):
    # no trips at all, or trips of less than 3km in total.
    sample_am = [seg for seg in segs if time_overlap(seg, tam_0, tam_1) and seg['type'] in ['trip', 'jump', 'jumptrip']]
    trips_at_least_3k_am = sum([seg_dist(seg) for seg in sample_am]) > 3.0
    flag_no_trip_am = not sample_am or not trips_at_least_3k_am

    # no trips at all, or trips of less than 3km in total.
    sample_pm = [seg for seg in segs if time_overlap(seg, tpm_0, tpm_1) and seg['type'] in ['trip', 'jump', 'jumptrip']]
    trips_at_least_3k_pm = sum([seg_dist(seg) for seg in sample_pm]) > 3.0
    flag_no_trip_pm = not sample_pm or not trips_at_least_3k_pm

    return flag_no_trip_am, flag_no_trip_pm


def flag_outstation(segs):

    # start day in Bangalore
    start_bangalore = False
    locs = [seg for seg in segs if seg['type'] == 'loc']
    if locs: start_bangalore = int(locs[0]['outstation_full']) == 0

    # end day in Bangalore
    end_bangalore = False
    locs = [seg for seg in segs if seg['type'] == 'loc']
    if locs: end_bangalore = int(locs[0]['outstation_full']) == 0

    time_out = sum([time_overlap(seg, 7, 22) for seg in segs if int(seg['outstation_full']) == 1])
    time_tot = sum([time_overlap(seg, 7, 22) for seg in segs])
    if time_tot == 0.0:
        return False, False
    elif time_out / time_tot > 0.5 and not (start_bangalore or end_bangalore):
        return True, True
    else:
        return False, False

    # flag_out_am = any([bool(int(seg['outstation_full'])) for seg in segs if time_overlap(seg, 7, 14)])
    # flag_out_pm = any([bool(int(seg['outstation_full'])) for seg in segs if time_overlap(seg, 14, 22)])
    # return flag_out_am, flag_out_pm


def compute_charge_dt(start_time, seg_dist, d0, max_rate_am):
    """
    ramp = 0 outside
    ramp = -1 first shoulder
    ramp = 2 peak
    ramp = 1 second shoulder
    ramp = 0 outside
    :param start_time:
    :param seg_dist:
    :param d0:
    :param charge_am:
    :return:
    """

    seg_dist = floor(10.0 * seg_dist) / 10.0

    start_time = ut.hf_hhmm(start_time)
    if d0 < start_time <= d0 + 1:
        ramp = -1
        rate = (start_time - d0) * max_rate_am
        rate = floor(2.0 * rate) / 2.0
        tot_charge = rate * seg_dist
        tot_charge = floor(tot_charge)

    elif d0 + 1 < start_time < d0 + 2:
        ramp = 2
        rate = max_rate_am
        tot_charge = rate * seg_dist
        tot_charge = floor(tot_charge)

    elif d0 + 2 <= start_time < d0 + 3:
        ramp = 1
        rate = (d0 + 3 - start_time) * max_rate_am
        rate = floor(2.0 * rate) / 2.0
        tot_charge = rate * seg_dist
        tot_charge = floor(tot_charge)
    else:
        ramp = 0
        rate = 0.0
        tot_charge = 0

    return ramp, rate, seg_dist, tot_charge


def compute_charges_dt_chain(respondent, mydate, buffer=1, treat_params=None, segs=None, chains=None):
    myuidp = respondent['uidp']

    # assert segs
    # assert structs
    print('a', end="")

    # load DT parameters
    tam0 = float(treat_params['DT AM start'])
    tpm0 = float(treat_params['DT PM start'])

    # high/low charge hardcoded as 24 and 12
    max_rate_am = float(treat_params['DT treat'])
    max_rate_pm = max_rate_am

    charge_max_daily = int(float(treat_params['DT max daily']))
    charge_no_data = int(float(treat_params['DT no data']))
    charge_no_trip = int(float(treat_params['DT no trip']))

    """ no data """
    flag_missing_data_am, flag_missing_data_pm, stats = flag_no_data(segs, tam0, tpm0)

    """ no trip - in +/-2h interval (7 hours) """
    flag_no_trip_am, flag_no_trip_pm = flag_no_trip(segs, tam0 - 2, tam0 + 5, tpm0 - 2, tpm0 + 5)

    """ outstation - using "full" definition = seg['d2b_p20'] > 18000 """
    flag_out_am, flag_out_pm = flag_outstation(segs)

    """ trips """
    trips_am = []
    trips_pm = []

    for chain in chains:
        if chain['trip'] == '1':
            chain['t01'] = chain['t_start_bg'].zfill(8)[:5]
            chain['t11'] = chain['t_stop_bg'].zfill(8)[:5]
            seg_dist_km = float(chain['plb']) / 1000

            # intervals for AM and PM trips
            am_min = tam0 - buffer
            am_max = min(tam0 + 3 + buffer, 14.5)
            pm_min = max(14.5, tpm0 - buffer)
            pm_max = tpm0 + 3 + buffer

            if time_overlap(chain, am_min, am_max):
                ramp, rate, my_seg_dist, tot_charge = compute_charge_dt(chain['t01'], seg_dist_km, tam0, max_rate_am)
                trips_am.append({
                    'ramp': ramp,
                    'rate': rate,
                    'charge': tot_charge,
                    'start': chain['t01'],
                    'duration': int(round(float(chain['dur_bg_mm']))),
                    'distance': my_seg_dist,
                    'id_new': chain['chain']
                })
            elif time_overlap(chain, pm_min, pm_max):
                ramp, rate, my_seg_dist, tot_charge = compute_charge_dt(chain['t01'], seg_dist_km, tpm0, max_rate_pm)
                trips_pm.append({
                    'ramp': ramp,
                    'rate': rate,
                    'charge': tot_charge,
                    'start': chain['t01'],
                    'duration': int(round(float(chain['dur_bg_mm']))),
                    'distance': my_seg_dist,
                    'id_new': chain['chain']
                })

    # be lenient: only mark flag_no_trip = TRUE if BOTH are true
    # (only do this if there are NO trips at all)
    if not (flag_no_trip_am and flag_no_trip_pm):
        flag_no_trip_am = False
        flag_no_trip_pm = False

    """ Compute charges """
    # AM
    if flag_missing_data_am:
        charge_am = charge_no_data
    elif flag_no_trip_am:
        charge_am = charge_no_trip
    else:
        charge_am = 0
        for trip in trips_am:
            charge_am += trip['charge']

    # PM
    if flag_missing_data_pm:
        charge_pm = charge_no_data
    elif flag_no_trip_pm:
        charge_pm = charge_no_trip
    else:
        charge_pm = 0
        for trip in trips_pm:
            charge_pm += trip['charge']

    # total daily
    charge_total_raw = charge_am + charge_pm
    charge_total = min(charge_max_daily, charge_am + charge_pm)

    if flag_out_am or flag_out_pm:
        charge_am = 0
        charge_pm = 0
        charge_total = 0

    results = {
        'uidp': myuidp,
        'date': mydate,
        'flag_no_data_am': flag_missing_data_am,
        'flag_no_data_pm': flag_missing_data_pm,
        'flag_no_trip_am': flag_no_trip_am,
        'flag_no_trip_pm': flag_no_trip_pm,
        'flag_out_am': flag_out_am,
        'flag_out_pm': flag_out_pm,
        'outstation': flag_out_am or flag_out_pm,
        'trips_am': trips_am,
        'trips_pm': trips_pm,
        'charge_am': charge_am,
        'charge_pm': charge_pm,
        'charge_total': charge_total,
        'charge_total_raw': charge_total_raw
    }
    return results


def compute_charge_area_chain(chain, segs_dict, points, area_center, area_radius):
    """
    more sophisticated version for analysis. Compute all intersections and if jump/normal
     -- for each segment, compute min distance
     -- keep all within or intersecting circle
     -- get minimum *intersecting* seg length
     -- get first time of intersection
    """

    chain['points'] = []
    # pprint(sorted([pt['idx'] for pt in points]))
    segs_list = eval(chain['segs'])
    for idx in segs_list:
        seg = segs_dict[idx]
        if seg['type'] in ['jump', 'jumptrip', 'trip', 'gap']:
            i_orig = int(seg['start'])
            i_dest = int(segs_dict[idx + 1]['start'])  # todo: check that this always exists

            # note: we select only NON-dropped points,
            # except if they are endpoints (it means they were low acc addeed back)
            chain['points'] += [{
                                   'time': pt['time'],  # pt['date'] + " " +
                                   'latitude': float(pt['latitude']),
                                   'longitude': float(pt['longitude'])
                               } for pt in points if i_orig <= int(pt['idx']) <= i_dest
                                and (pt['drop'] == '0' or i_orig == int(pt['idx']) or int(pt['idx']) == i_dest)]

    if len(chain['points']) <= 1:
        print(chain['uidp'])
        print(chain['date'])
        print(chain['chain'])
        print(chain['points'])
        raise ArithmeticError

    # for each edge in the trip, compute the minimum distance from the area center to the edge
    idx_within = []
    for idx, pt in enumerate(chain['points'][:-1]):
        pt_next = chain['points'][idx + 1]
        start = (pt['latitude'], pt['longitude'])
        end = (pt_next['latitude'], pt_next['longitude'])
        pt['dist2edge'], _ = utg.pnt2line(pnt=area_center, start=start, end=end)
        pt['edge_len'] = great_circle(start, end).meters
        if pt['dist2edge'] < area_radius:
            idx_within.append(idx)

    min_dist = min([pt['dist2edge'] for pt in chain['points'][:-1]])
    intersect = min_dist < area_radius
    min_edge = -1
    if idx_within:
        min_edge = min([chain['points'][idx]['edge_len'] for idx in idx_within
                        if chain['points'][idx]['edge_len'] > 0.0])

    first_time = ''
    if idx_within:
        idx = idx_within[0]
        first_time = chain['points'][idx]['time']
    return intersect, min_dist, first_time, min_edge


def compute_charges_area_chain(ROOT_PATH, respondent, mydate, chain_th=15, treat_params=None):
    myuidp = respondent['uidp']

    # mydate YYYY-MM-DD

    # load chain data
    chains_data_file = ROOT_PATH + 'data/coded_gps_byid/' + myuidp + '/segs_chain_' + str(chain_th) + '.csv'
    chains = ut.csv2dict(chains_data_file)
    chains = [chain for chain in chains if chain['date'] == mydate]

    # load (meta) trip data
    trip_data_file = ROOT_PATH + 'data/coded_gps_byid/' + myuidp + '/trips.csv'
    segs = ut.csv2dict(trip_data_file)
    segs = [seg for seg in segs if seg['date'] == mydate and seg['type'] != 'end']
    segs_dict = {int(seg['id']): seg for seg in segs}

    # todo: achieve differently?
    for idx, seg in enumerate(segs):
        if seg['type'] in ['trip', 'jump', 'jumptrip']:
            if idx + 1 < len(segs):
                seg['end'] = segs[idx + 1]['start']
            else:
                seg['end'] = 2000000
    # pprint_segs(segs)

    # load points data - points_L0202093440_2017-02-18.csv
    points_data_file = ROOT_PATH + 'data/coded_gps_byid/' + myuidp + '/points/points_' + myuidp + '_' + mydate + '.csv'
    my_points = ut.csv2dict(points_data_file)
    my_points = [point for point in my_points]  # if point['drop'] == '0'

    # check that AREA treatment exists
    if respondent['atreat'] == '0 No treatment':
        return False

    # load AREA parameters
    area_center = (float(treat_params['alat']),
                       float(treat_params['alng']),)
    area_radius = int(treat_params['aradius'])
    area_charge = int(treat_params['acharge'])

    """ no data """
    flag_missing_data_am, flag_missing_data_pm, _ = flag_no_data(segs, tam0=8.5, tpm0=17.5)

    """ no trip """
    flag_no_trip_am, flag_no_trip_pm = flag_no_trip(segs, 7, 14, 14, 21)

    """ outstation - using "full" definition = seg['d2b_p20'] > 18000 """
    flag_out_am, flag_out_pm = flag_outstation(segs)

    # be lenient: only mark flag_no_trip = TRUE if BOTH are true
    # (only do this if there are NO trips at all)
    if not (flag_no_trip_am and flag_no_trip_pm):
        flag_no_trip_am = False
        flag_no_trip_pm = False

    """ Check if any of the trips intersects the area """
    intersect_trips = []
    # morning chains
    for chain in chains:
        if chain['trip'] == '1':
            # save t01 and t11
            chain['t01'] = chain['t_start_bg'].zfill(8)[:5]
            chain['t11'] = chain['t_stop_bg'].zfill(8)[:5]

            # intersection
            intersect, min_dist, time_intersect, min_edge = \
                compute_charge_area_chain(chain=chain, segs_dict=segs_dict, points=my_points,
                                          area_center=area_center, area_radius=area_radius)
            charge = intersect * area_charge
            # save
            intersect_trips.append({
                'id_new': chain['chain'],
                'intersect': int(intersect),
                'time_intersect': time_intersect,
                'min_dist': min_dist,
                'min_edge': min_edge,
                'charge': charge
            })

            # save to chain!
            chain['intersect'] = intersect

    # any AM/PM intersection?
    intersects_am = [chain for chain in chains if chain['trip'] == '1' and time_overlap(chain, 7, 14) > 0.0 and chain['intersect']]
    intersects_pm = [chain for chain in chains if chain['trip'] == '1' and time_overlap(chain, 14, 21) > 0.0 and chain['intersect']]
    if intersects_am: intersect_am = {'intersect': True, 'charge': area_charge}
    else:             intersect_am = {'intersect': False, 'charge': 0}
    if intersects_pm: intersect_pm = {'intersect': True, 'charge': area_charge}
    else:             intersect_pm = {'intersect': False, 'charge': 0}

    """ Compute charges """
    # AM
    if flag_missing_data_am:
        charge_am = int(respondent['anodata'])
    elif flag_no_trip_am:
        charge_am = int(respondent['anotrip'])
    else:
        charge_am = intersect_am['charge']

    # PM
    if flag_missing_data_pm:
        charge_pm = int(respondent['anodata'])
    elif flag_no_trip_pm:
        charge_pm = int(respondent['anotrip'])
    else:
        charge_pm = intersect_pm['charge']

    # total daily
    charge_total = min(int(respondent['amaxdaily']), charge_am + charge_pm)

    if flag_out_am or flag_out_pm:
        charge_am = 0
        charge_pm = 0
        charge_total = 0

    """ Generate daily summary """
    response = {
        'uidp': myuidp,
        'date': mydate,
        'flag_no_data_am': flag_missing_data_am,
        'flag_no_data_pm': flag_missing_data_pm,
        'flag_no_trip_am': flag_no_trip_am,
        'flag_no_trip_pm': flag_no_trip_pm,
        'flag_out_am': flag_out_am,
        'flag_out_pm': flag_out_pm,
        'outstation': flag_out_am or flag_out_pm,
        'intersect_trips': intersect_trips,
        'charge_am': charge_am,
        'charge_pm': charge_pm,
        'intersect_am': int(intersect_am['intersect']),
        'intersect_pm': int(intersect_pm['intersect']),
        'charge_intersect_am': int(intersect_am['charge']),
        'charge_intersect_pm': int(intersect_pm['charge']),
        'charge_total': charge_total
    }

    return response


def compute_charges(ROOT_PATH, chain_th=15):
    """ This script: compute charges for everyone in the analysis sample, for all days. Store. """

    """ paths # read sample """
    treatment_roster_file = ROOT_PATH + 'data/treatment/treatment roster noPII.csv'
    treat_roster = ut.csv2dict(treatment_roster_file)
    treat_roster = [line for line in treat_roster if line['meeting'] == 'done']
    # treat_roster = treat_roster[8:12]  # testing
    # print(treat_roster)

    all_charges = []  # daily
    all_trips_dt = []
    all_trips_area = []
    for resp_idx, resp in enumerate(treat_roster):
        uidp = resp['uidp']
        print("")
        print(str(resp_idx+1) + "/" + str(len(treat_roster) + 1) + " Processing uidp " + uidp)

        """ read chains """
        chains_file = ROOT_PATH + 'data/coded_gps_byid/' + uidp + '/segs_chain_' + str(chain_th) + '.csv'
        all_chains = ut.csv2dict(chains_file)

        """ read simple trips """
        trips_file = ROOT_PATH + 'data/coded_gps_byid/' + uidp + '/trips.csv'
        all_segs = ut.csv2dict(trips_file)
        for seg in all_segs:
            seg['path_length_inbang_tot'] = seg['path_length_inbang']

        """ all dates for this uidp """
        all_dates = set([datetime.datetime.strptime(line['date'], '%Y-%m-%d').date() for line in all_chains])
        all_dates = sorted(list(all_dates))

        if len(all_dates) == 0:
            continue

        # add all missing intermediate dates, except weekends
        min_date = min(all_dates)
        max_date = max(all_dates)
        all_dates = [min_date + datetime.timedelta(days=i) for i in range((max_date - min_date).days + 1)]
        # all_dates = [str(mydate) for mydate in all_dates if mydate.weekday() not in [5, 6]]  # 0 = Monday
        all_dates = [str(mydate) for mydate in all_dates]  # 0 = Monday
        assert str(max_date) in all_dates or max_date.weekday() in [5, 6]

        """ read treat_params from treat prep file """
        dt_treat_params_file = ROOT_PATH + 'data/treatment/treatment_parameters/' + uidp + '/dt_treat_parameters.csv'
        dt_treat_params = ut.csv2dict(dt_treat_params_file)

        """ Loop over all dates in data """
        for mydate in all_dates:

            """ DT Charges """
            # Loop over scenarios (low/high balance target, low/high charges)
            for i in [0, 1, 2, 3]:
                my_dt_params = dt_treat_params[i]
                bal_target = my_dt_params['bal_target_name']
                h = my_dt_params['h_name']

                treat_params_dt = {
                    'DT AM start': my_dt_params['tam_d0'],
                    'DT PM start': my_dt_params['tpm_d0'],
                    'DT treat': my_dt_params['h'],
                    # 'DT treat': 24,
                    'DT max daily': my_dt_params['daily_charge_maximum'],
                    'DT no data': my_dt_params['charge_nodata'],
                    'DT no trip': my_dt_params['charge_notrip']  # todo: why was this switched to zero?
                }

                """ Compute DT Charges at DAILY AND TRIP level """
                segs = [line for line in all_segs if line['date'] == str(mydate) and line['type'] != 'end']
                chains = [line for line in all_chains if line['date'] == str(mydate)]
                charges_dt = compute_charges_dt_chain(resp, mydate, buffer=4,
                                                      treat_params=treat_params_dt,
                                                      segs=segs, chains=chains)

                """ Save DAILY LEVEL """
                charges_dt = ut.bool2int(charges_dt)
                charges_dt['type'] = 'dt'
                charges_dt['bal_target'] = bal_target
                charges_dt['h'] = h
                charges_dt['charge_nodata'] = my_dt_params['charge_nodata']
                charges_dt['charge_notrip'] = my_dt_params['charge_notrip']
                charges_dt['charge_maxday'] = my_dt_params['daily_charge_maximum']

                # save charges, save trips
                all_charges.append(charges_dt)

                """ Save TRIP LEVEL """
                # only save trips once
                if h == 'high' and bal_target == 'high':
                    trips_am = charges_dt.pop('trips_am', [])
                    trips_pm = charges_dt.pop('trips_pm', [])

                    for trip in trips_am + trips_pm:
                        trip['uidp'] = uidp
                        trip['date'] = mydate

                    # save trips
                    all_trips_dt += trips_am + trips_pm

            """ AREA Charges """
            # pprint(resp)

            # if resp['A treat'] != '0 No treatment':
            if resp['atreat'] != '0 No treatment':
                # a single scenario for everyone (these are hypothetical charges)
                treat_params_area = {
                    'alat': resp['alat'],
                    'alng': resp['alng'],
                    'aradius': resp['aradius'],
                    'acharge': 160
                }

                """ Compute AREA Charges at DAILY AND TRIP level """
                charges_area_new = compute_charges_area_chain(ROOT_PATH=ROOT_PATH, respondent=resp,
                                                              mydate=mydate, chain_th=chain_th,
                                                              treat_params=treat_params_area)

                """ Save DAILY LEVEL """
                charges_area_new = ut.bool2int(charges_area_new)
                charges_area_new['type'] = 'area_new'
                all_charges.append(charges_area_new)

                """ Save TRIP LEVEL """
                intersect_trips = charges_area_new.pop('intersect_trips', [])
                for trip in intersect_trips:
                    trip['uidp'] = uidp
                    trip['date'] = mydate
                all_trips_area += intersect_trips

    """ Save to file """
    print("\nfinished processing charges")

    """ DT - TRIPS """
    all_trips_dt_file = ROOT_PATH + "data/coded_gps/charges_trips_dt_" + str(chain_th) + ".csv"
    ut.dictlist2csv(all_trips_dt, all_trips_dt_file)
    print("Saved DT trips")

    """ AREA - TRIPS """
    all_trips_area_file = ROOT_PATH + "data/coded_gps/charges_trips_area_" + str(chain_th) + ".csv"
    ut.dictlist2csv(all_trips_area, all_trips_area_file)
    print("Saved AREA trips")

    """ DT and AREA - DAILY"""
    all_charges_file = ROOT_PATH + "data/coded_gps/charges_daily_all_" + str(chain_th) + ".csv"
    print("Sorting... ")
    all_charges = sorted(all_charges, key=itemgetter('uidp', 'date', 'type'))
    ut.dictlist2csv(all_charges, all_charges_file)
    print("Saved charges")


if __name__ == '__main__':
    ROOT_PATH = 'C:/bang_cp_paper/'
    ROOT_PATH = 'C:/Users/Gabriel Kreindler/Dropbox/projects/bang_cp_paper/replication_package/'

    t0 = time.time()
    compute_charges(ROOT_PATH=ROOT_PATH, chain_th=15)
    compute_charges(ROOT_PATH=ROOT_PATH, chain_th=30)
    compute_charges(ROOT_PATH=ROOT_PATH, chain_th=60)
    t1 = time.time()
    print("Time for segments coding (all subjects, all dates): " + str(floor(t1 - t0)))
