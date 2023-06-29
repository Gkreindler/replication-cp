import datetime
import os
import time
import numpy
import detect
from operator import itemgetter
import pointclass as pc
import tripcoding_functions as tp_fn
import analysis_prepare_sample as aps
import utils as ut
from joblib import Parallel, delayed
import multiprocessing
import math
import gps_data_workflow as gpsdw
import clean_conflicting_copies_dropbox as cccd

numpy.set_printoptions(precision=3, suppress=True)


def deviceid_coding(ROOT_PATH):

    # read from files by date, minimal processing, dump by deviceid
    gpsdw.read_dump_deviceid(ROOT_PATH=ROOT_PATH, line_buffer_size=2000000)

    # remove duplicates and sort by time
    gpsdw.clean_sort_deviceid(ROOT_PATH=ROOT_PATH)

    """ refresh deviceid-devicename crosswalk """
    # clean dropbox duplicates of deviceid's;j
    cccd.clean_subfolder(ROOT_PATH + 'data/raw_gps/', readonly=False)
    cccd.clean_tree(ROOT_PATH + 'data/', readonly=False)

    # check for any deviceids that appear in CTO but not in the data (nor online after ping)
    gpsdw.crosswalk_check(ROOT_PATH=ROOT_PATH)

    return True


def tp_points_onedate(points_list, date, uidp, deviceid,
                      drop_color=False, debug=False,
                      pts2inspect=None, OUTPATHSTUB=None, FOLDEROUT_DEBUG=None,
                      result="", prefix=""):
    out_text = prefix + uidp + ' ' + str(date) + ' n points ' + str(len(points_list))
    print(prefix + out_text)
    result += out_text + '\n'

    if not pts2inspect:
        pts2inspect = []

    """ high radius """
    detect.detect_highradius(points_list, min_accuracy_radius=400)
    print(prefix + "POINT CODE: high radius drop - done")

    """ refresh angle, dist2prev/next and times """
    pc.point_stats(points_list, debug=debug)
    print(prefix + "POINT CODE: stats - done")

    """ remove crowded outliers """
    detect.detect_crowded_outliers(points_list)
    print(prefix + "POINT CODE: detected crowded outliers - done")

    """ sharp angles """
    for _ in range(2):
        # refresh angle, dist2prev/next and times
        pc.point_stats(points_list)
        print(prefix + "POINT CODE: stats - done")

        # sharp angles
        detect.detect_sharp_angle(points_list)
        print(prefix + "POINT CODE: sharp jump drop - done")

    """ DETECT ANGLE with multiple points """
    for ir in range(2):
        print(prefix + "---- sharp mult angle ---- " + str(ir))
        pts2inspect += detect.detect_sharp_multiple(points_list, go_fwd=True, debug=debug, prefix=prefix)
        print(prefix + "---- sharp mult angle in reverse ---- " + str(ir))
        pts2inspect += detect.detect_sharp_multiple(points_list, go_fwd=False, debug=debug, prefix=prefix)  # in reverse

    for ir in range(2):
        """ DETECT LAZY """
        print(prefix + "---- lazy ---- " + str(ir))
        pts2inspect += detect.detect_lazy(points_list, go_fwd=True)
        print(prefix + "---- lazy in reverse ---- " + str(ir))
        pts2inspect += detect.detect_lazy(points_list, go_fwd=False)  # in reverse

    for ir in range(2):
        """ DROP 1st and last point if it's high speed """
        print(prefix + "---- first ---- " + str(ir))
        pts2inspect += detect.drop_first_high_speed(points_list, go_fwd=True)
        print(prefix + "---- last ---- " + str(ir))
        pts2inspect += detect.drop_first_high_speed(points_list, go_fwd=False)

    """ one more sharp angles """
    # refresh angle, dist2prev/next and times
    pc.point_stats(points_list)
    print(prefix + "POINT CODE: stats - done")
    # sharp angles
    detect.detect_sharp_angle(points_list)
    print(prefix + "POINT CODE: sharp jump drop - done")

    """ DETECT OUTSTANDING LARGE SPEED """
    print(prefix + "---- high speed ----")
    pts2inspect += detect.detect_high_speed(points_list, go_fwd=True)
    print(prefix + "---- high speed in reverse ----")
    pts2inspect += detect.detect_high_speed(points_list, go_fwd=False)  # in reverse

    """ TRIP CODING """
    # add score_after
    detect.location_score(points_list)
    print(prefix + "POINT CODE: location score - done")

    detect.detect_trips(points_list)
    print(prefix + "POINT CODE: detect trips - done")

    # refresh angle, dist2prev/next and times
    pc.point_stats(points_list)
    print(prefix + "POINT CODE: stats - done")

    """ WRITE TO FILE """
    if points_list:
        print(prefix + "POINT CODE: writing to file")

        if not FOLDEROUT_DEBUG:
            # prepare to write output - check that relevant folders exists
            FOLDEROUT = OUTPATHSTUB + uidp
            if not os.path.exists(FOLDEROUT):
                try:
                    os.makedirs(FOLDEROUT)
                except:
                    pass
            if not os.path.exists(FOLDEROUT + "/points/"):
                try:
                    os.makedirs(FOLDEROUT + "/points/")
                except:
                    pass
            PATHOUT = FOLDEROUT + "/points/points_" + uidp + "_" + date.strftime("%Y-%m-%d") + ".csv"
        else:
            # write pts2inspect
            # pts2inspect_file = FOLDEROUT_DEBUG + "/points2inspect.csv"
            # ut.dictlist2csv(pts2inspect, pts2inspect_file)
            PATHOUT = FOLDEROUT_DEBUG + "/points/points_" + uidp + "_" + date.strftime("%Y-%m-%d") + ".csv"

        # write to file
        ut.pointlist2csv(points_list, PATHOUT, drop_color=drop_color)
        out_text = prefix + "POINT CODE: Points file successfully written to disk: " + \
                   deviceid + ", " + str(date) + " UIDP " + uidp + '\n'
        print(out_text)
        result += out_text
    else:
        out_text = prefix + "POINT CODE: Empty list of points: " + \
                   deviceid + ", " + str(date) + " - no file created " + " UIDP " + uidp + '\n'
        print(out_text)
        result += out_text

    return result, pts2inspect


def tp_points_parallel(uidp_date_dict,
                       drop_color=True, debug=False, num_cores=4, n_points_threshold=300000,
                       FOLDEROUT_DEBUG=None, ROOT_PATH=None, OUTPATHSTUB=None):
    """ """

    """ Clean, detect trips and write to file """
    if FOLDEROUT_DEBUG:
        pts2inspect_file = FOLDEROUT_DEBUG + "/points2inspect.csv"
        pts2inspect = ut.csv2dict(pts2inspect_file)
    else:
        pts2inspect = []

    """ Load points data for specific users """
    print("POINT CODE: read all BTR data", end="")
    allpoints = {}
    uidp_date_list = []
    n_points = 0
    all_results = []
    prefix_stage = 1
    skipped_devids = set()
    skipped_errors = []
    for idx_uidp, uidp in enumerate(sorted(uidp_date_dict.keys())):
        deviceid = uidp_date_dict[uidp]['deviceid']
        dates = sorted(list(uidp_date_dict[uidp]['dates']))
        prefix = uidp_date_dict[uidp]['prefix']

        PATHIN = ROOT_PATH + "data/raw_gps/coded_btr_" + deviceid + ".csv"
        try:
            allpoints[uidp] = list(pc.generate_points(PATHIN))
            print("READ " + uidp)

            # drop points without latitude or longitude (i.e. equal to zero)
            # use abs() > 0.01 as tolerance
            allpoints[uidp] = [point for point in allpoints[uidp] if abs(point.latitude) > 0.01 and abs(point.longitude) > 0.01]
            print("dropped points without lat and long - done")
            n_points += len(allpoints[uidp])

            # split by day, and sort within day by timestamp
            allpoints[uidp] = pc.sort_by_date(allpoints[uidp], dates=dates)
            print(prefix + "POINT CODE: sort by date - done")

            # create list of uidp-date pairs
            dates = [date for date in dates if date in allpoints[uidp]]
            for idx_date, date in enumerate(dates):
                uidp_date_list.append({
                    'uidp': uidp,
                    'deviceid': deviceid,
                    'date': date,
                    'prefix': uidp_date_dict[uidp]['prefix'] +
                              str(idx_uidp + 1) + "/" + str(len(uidp_date_dict)) + ", " +
                              str(idx_date + 1) + "/" + str(len(dates)) + " "
                })

            """ Run trip-coding every N_POINTS =  """
            if n_points >= n_points_threshold or idx_uidp == len(uidp_date_dict.keys()) - 1:
                print("POINT CODE: clean, detect trips for " + str(n_points) + " number of points")

                prefix_stage_str = 'stage ' + str(prefix_stage) + ' '

                if allpoints:
                    # parallel computing with virtual cores
                    # num_cores = multiprocessing.cpu_count()
                    if num_cores > 1:
                        print("Running PARALLEL with " + str(num_cores) + " cores")
                        parallel_output = Parallel(n_jobs=num_cores) \
                            (delayed(tp_points_onedate) \
                                 (points_list=allpoints[item['uidp']][item['date']],
                                  date=item['date'], uidp=item['uidp'], deviceid=item['deviceid'],
                                  drop_color=drop_color, debug=debug, pts2inspect=pts2inspect,
                                  OUTPATHSTUB=OUTPATHSTUB, FOLDEROUT_DEBUG=FOLDEROUT_DEBUG,
                                  prefix=prefix_stage_str + item['prefix']) \
                             for item in uidp_date_list)

                        all_results += [output[0] for output in parallel_output]
                        pts2inspect += [item for output in parallel_output
                                        for item in output[1]]
                    else:
                        print("Running NON Parallel")
                        for item in uidp_date_list:
                            uidp = item['uidp']
                            date = item['date']
                            deviceid = item['deviceid']
                            result_output, pts2inspect_output = \
                                tp_points_onedate(points_list=allpoints[uidp][date],
                                                  date=date, uidp=uidp, deviceid=deviceid,
                                                  drop_color=drop_color, debug=debug, pts2inspect=pts2inspect,
                                                  OUTPATHSTUB=OUTPATHSTUB, FOLDEROUT_DEBUG=FOLDEROUT_DEBUG,
                                                  prefix=prefix_stage_str + item['prefix'])
                            all_results += result_output
                            pts2inspect += pts2inspect_output

                """ Reset the list """
                n_points = 0
                uidp_date_list = []
                allpoints = {}
                prefix_stage += 1
        except Exception as e:
            print("Problem with data download " + deviceid)
            skipped_devids.add(deviceid)
            skipped_errors.append({
                'deviceid': deviceid,
                'error': str(e)
            })

    """ write pts2inspect """
    if FOLDEROUT_DEBUG:
        pts2inspect_file = FOLDEROUT_DEBUG + "/points2inspect.csv"
        ut.dictlist2csv(pts2inspect, pts2inspect_file)
    return " ".join(all_results)


def run_points_coding(ROOT_PATH, dates, uidp_list=None, OUTPATHSTUB="data/coded_gps_byid/",
                      num_cores=4, n_points_threshold=300000, skip_existing_folders=False):

    PATH_UIDP_DEVICEID = ROOT_PATH + "data/coded_cto/crosswalk_uidp_deviceid.csv"
    crosswalk_list = ut.csv2dict(PATH_UIDP_DEVICEID)

    if uidp_list:
        crosswalk_list = [line for line in crosswalk_list if line['uidp'] in uidp_list]

    uidp_date_dict = {}
    for idx, line in enumerate(crosswalk_list):
        try:
            uidp = line['uidp']
            deviceid = line['deviceid']

            skip_this_uidp = False
            if skip_existing_folders:
                uidp_path = ROOT_PATH + OUTPATHSTUB + uidp
                skip_this_uidp = os.path.isdir(uidp_path)

            if not skip_this_uidp:
                assert uidp not in uidp_date_dict
                uidp_date_dict[uidp] = {
                    'uidp': uidp,
                    'dates': {date for date in dates},
                    'deviceid': deviceid,
                    'prefix': ''
                }
            else:
                print("skipping " + uidp + " because the folder already exists")
        except Exception as e:
            raise e
    print("POINTS PREP DONE")

    """ run """
    result = []
    result += tp_points_parallel(uidp_date_dict,
                                ROOT_PATH=ROOT_PATH, OUTPATHSTUB=ROOT_PATH + OUTPATHSTUB,
                                num_cores=num_cores, n_points_threshold=n_points_threshold)

    return result


def run_seg_coding(ROOT_PATH, uidps,
                   chains_max_loc_dur_mm=30,
                   PATHIN_STUB=None,
                   PATHOUT_STUB=None,
                   debug=False, silent=False, num_cores=1):
    """
    by uidp, read all points (by date), organize into segments, and then improve segments
    :param uidps:
    :return:
    """

    if not PATHIN_STUB:
        PATHIN_STUB = ROOT_PATH + 'data/raw_gps/'
    if not PATHOUT_STUB:
        PATHOUT_STUB = ROOT_PATH + 'data/coded_gps_byid/'

    PATH_UIDP_DEVICEID = ROOT_PATH + "data/coded_cto/crosswalk_uidp_deviceid.csv"
    crosswalk_dict = ut.csv2dict(PATH_UIDP_DEVICEID)
    if uidps:
        crosswalk_dict = [subject for subject in crosswalk_dict if subject['uidp'] in uidps]

    t0 = time.time()

    # parallel
    if num_cores > 1:
        num_cores = min(num_cores, multiprocessing.cpu_count())
        print("Number of virtual cores to use: " + str(num_cores))
        Parallel(n_jobs=num_cores) \
            (delayed(tp_fn.tp_segments) \
                 (deviceid=subject["deviceid"],
                  uidp=subject["uidp"],
                  name=subject.get("name", "no_name"),
                  chains_max_loc_dur_mm=chains_max_loc_dur_mm,
                  PATHIN_STUB=PATHIN_STUB,
                  PATHOUT_STUB=PATHOUT_STUB,
                  ROOT_PATH=ROOT_PATH,
                  debug=debug,
                  silent=silent,
                  prefix=subject["uidp"] + " " + str(idx+1) + '/' + str(len(crosswalk_dict)) + ". ") \
             for idx, subject in enumerate(crosswalk_dict))

    else:
        # sequential
        for idx, subject in enumerate(crosswalk_dict):
            deviceid, uidp, name = subject["deviceid"], subject["uidp"], subject.get("name", "no_name")
            print("Processing segments " + str([deviceid, uidp, name]))
            tp_fn.tp_segments(deviceid=deviceid, uidp=uidp, name=name,
                              chains_max_loc_dur_mm=chains_max_loc_dur_mm,
                              PATHIN_STUB=PATHIN_STUB,
                              PATHOUT_STUB=PATHOUT_STUB,
                              ROOT_PATH=ROOT_PATH,
                              debug=debug,
                              prefix=subject["uidp"] + " " + str(idx + 1) + '/' +
                                     str(len(crosswalk_dict)) + ". ")

    t1 = time.time()
    print("Time for segments coding (all subjects, all dates): " + str(math.floor(t1 - t0)))


def analysis_segment_coding(ROOT_PATH, uidps=None, silent=False, n_cores=8):
    PATHIN_STUB = ROOT_PATH + 'data/coded_gps_byid/'
    PATHOUT_STUB = ROOT_PATH + 'data/coded_gps_byid/'
    debug = False

    chains_max_loc_dur_mm = 30

    # EXPERIMENT uidps
    TREAT_PATH = ROOT_PATH + 'treatment/treatment roster noPII.csv'
    treat_roster = ut.csv2dict(TREAT_PATH)
    if not uidps:
        uidps = [line['uidp'] for line in treat_roster if line['meeting'] == 'done']
    # uidps = ['L0207085326']
    # debug = True

    num_cores = min(n_cores, len(uidps))
    run_seg_coding(ROOT_PATH=ROOT_PATH, uidps=uidps,
                   chains_max_loc_dur_mm=chains_max_loc_dur_mm,
                   PATHIN_STUB=PATHIN_STUB,
                   PATHOUT_STUB=PATHOUT_STUB,
                   debug=debug, num_cores=num_cores, silent=silent)


def quality_measures(ROOT_PATH, uidp_list):
    """ Compute a general proxy for data quality each day"""

    quality_stats = []

    # all dates
    first_day = datetime.date(year=2017, month=2, day=7)
    last_day = datetime.date(year=2017, month=9, day=11)
    all_dates = [first_day + datetime.timedelta(days=i) for i in range((last_day - first_day).days + 1)]

    for idx, uidp in enumerate(uidp_list):
        print(str(idx + 1) + "/" + str(len(uidp_list)) + "Processing " + uidp)

        trips_file = ROOT_PATH + "data/coded_gps_byid/" + uidp + "/trips.csv"
        segs = ut.csv2dict(trips_file)
        segs = [seg for seg in segs if seg['type'] != 'end']

        # all dates after recruitment
        mydates = [date for date in all_dates if date >= datetime.date(year=2017, month=int(uidp[1:3]), day=int(uidp[3:5]))]

        for mydate in mydates:
            mydate_str = mydate.strftime("%Y-%m-%d")
            segs_today = [seg for seg in segs if seg['date'] == mydate_str]
            qm_dict = tp_fn.compute_qm(uidp=uidp, segs_today=segs_today, mydate=mydate_str)
            quality_stats.append(qm_dict)

    quality_stats = sorted(quality_stats, key=(itemgetter('uidp', 'date')))
    ut.dictlist2csv(dictlist=quality_stats, pathout=ROOT_PATH + 'data/coded_gps/quality_all.csv')

    return True


def count_n_raw_points():
    # after the data is coded, count the total number of GPS points. Just for fun. 22,434,466
    # result:

    ALL_UIDPS_PATH = ROOT_PATH + 'data/coded_cto/crosswalk_uidp_deviceid.csv'
    uidp_list = ut.csv2dict(ALL_UIDPS_PATH)
    print(len(uidp_list))
    uidp_list = [line['uidp'] for line in uidp_list if
                 os.path.isfile(ROOT_PATH + 'data/raw_gps/coded_btr_' + line['deviceid'] + '.csv')]
    # uidp_list = uidp_list[0:2]
    print(len(uidp_list))
    np_total = 0
    for idx, uidp in enumerate(uidp_list):
        print(idx, uidp)
        file = ROOT_PATH + 'data/coded_gps_byid/' + uidp + '/trips.csv'
        temp = ut.csv2dict(file)
        np_total += sum([int(line['np']) for line in temp if line['np']])
    print(np_total)
    exit(0)

if __name__ == '__main__':
    # ROOT_PATH = 'C:/Users/Gabriel Kreindler/Dropbox/projects/bang_cp_paper/replication_package/'
    ROOT_PATH = '<YOUR_PATH_TO_REPLICATION_PACKAGE>'
    PATH_STUB = ROOT_PATH + 'data/raw_gps/'  # raw GPS data

    """ CODE RAW DATA BY DEVICEID """
    deviceid_coding(ROOT_PATH=ROOT_PATH)

    """ LOAD ALL UIDPS """
    ALL_UIDPS_PATH = ROOT_PATH + 'data/coded_cto/crosswalk_uidp_deviceid.csv'
    uidp_list = ut.csv2dict(ALL_UIDPS_PATH)
    print(len(uidp_list))
    uidp_list = [line['uidp'] for line in uidp_list if os.path.isfile(ROOT_PATH + 'data/raw_gps/coded_btr_' + line['deviceid'] + '.csv')]
    # uidp_list = uidp_list[0:2]
    print(len(uidp_list))
    print(uidp_list)

    # for processing one UID at a time
    # uidp_list = ['L0419164627']

    """ POINTS """
    skip_existing_folders = False
    d1 = datetime.date(2017, month=2, day=7)
    d2 = datetime.date(2017, month=9, day=18)  # inclusive
    date_list = [d1 + datetime.timedelta(days=i) for i in range((d2 - d1).days + 1)]
    run_points_coding(ROOT_PATH=ROOT_PATH, dates=date_list, OUTPATHSTUB="data/coded_gps_byid/", uidp_list=uidp_list,
                      num_cores=24, n_points_threshold=300000, skip_existing_folders=skip_existing_folders)

    """ SEGMENTS """
    n_hundreds = math.ceil(len(uidp_list)/100) + 1
    print(n_hundreds)
    for i in range(n_hundreds):
        print((100 * i))
        print((100 * (i+1)))
        uidp_list_chunk = uidp_list[(100 * i):(100 * (i+1))]
        if uidp_list_chunk:
            analysis_segment_coding(ROOT_PATH=ROOT_PATH, uidps=uidp_list_chunk, silent=True, n_cores=24)

    """ CHAINS """
    PATHIN = ROOT_PATH + "data/coded_gps_byid/"
    for idx, uidp in enumerate(uidp_list):
        print(str(idx) + '/' + str(len(uidp_list)))
        tp_fn.tp_chains(ROOT_PATH=ROOT_PATH, uidp=uidp, max_loc_dur_mm=60, PATHIN=PATHIN, debug=False)
        tp_fn.tp_chains(ROOT_PATH=ROOT_PATH, uidp=uidp, max_loc_dur_mm=30, PATHIN=PATHIN, debug=False)
        tp_fn.tp_chains(ROOT_PATH=ROOT_PATH, uidp=uidp, max_loc_dur_mm=15, PATHIN=PATHIN, debug=False)

    """ UPDATE CHAIN IDs """
    aps.update_chain_in_struct(uidps=uidp_list, PATHIN_STUB=PATHIN, chain_th=60)
    aps.update_chain_in_struct(uidps=uidp_list, PATHIN_STUB=PATHIN, chain_th=30)
    aps.update_chain_in_struct(uidps=uidp_list, PATHIN_STUB=PATHIN, chain_th=15)

    """ PREPARE SIMPLE ANALYSIS SAMPLE """
    aps.prepare_sample_complete(ROOT_PATH=ROOT_PATH, chain_th=15)
    aps.prepare_sample_complete(ROOT_PATH=ROOT_PATH, chain_th=30)
    aps.prepare_sample_complete(ROOT_PATH=ROOT_PATH, chain_th=60)

    aps.prepare_sample_experiment(ROOT_PATH=ROOT_PATH, chain_th=15)
    aps.prepare_sample_experiment(ROOT_PATH=ROOT_PATH, chain_th=30)
    aps.prepare_sample_experiment(ROOT_PATH=ROOT_PATH, chain_th=60)

    """ Quality Measures """
    quality_measures(ROOT_PATH=ROOT_PATH, uidp_list=uidp_list)

