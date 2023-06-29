import utils as ut
import codecs
import re
import sys
import gc
import glob
from datetime import datetime
from datetime import timedelta
from operator import itemgetter


def clean_apostrophe_file(file_path, replace=False):
    with codecs.open(file_path, "r", encoding='utf-8-sig') as file:
        file_data = file.read()
    # print(file_data)

    file_data_new = re.sub(r"([^\n,])\"([^\n,])", r"\1 \2", file_data)
    # file_data_new = re.sub(r"\",([^\n\"])", r" \1", file_data_new)
    # file_data_new = re.sub(r"([^\n\"]),\"", r"\1 ", file_data_new)

    if replace:
        with codecs.open(file_path, "w", encoding='utf-8-sig') as file:
            file.write(file_data_new)
    else:
        print(file_data_new)


def code_btr_data(lines, uidps, names, deviceid_replace, debug=False):
    """
    Do the basic coding and cleaning
    :param lines:
    :param uidps:
    :param names:
    :return:
    """

    time_format = '%a %b %d %Y %H:%M:%S'
    print("processing " + str(len(lines)) + " lines of data.")

    # line by line
    i = 0
    ids_to_delete = []
    for my_id in lines:

        # logging
        i += 1
        if i % 100000 == 1 and debug:
            print("<>", end="", flush=True)

        # latitude and longitude
        assert (lines[my_id]['latitude'] == '0' and lines[my_id]['longitude'] == '0') or \
               (lines[my_id]['latitude'] != '0' and lines[my_id]['longitude'] != '0')

        # change NEW deviceids to OLD ones (when a respondent changes his/her phone, etc)
        devid_temp = lines[my_id]['deviceid']
        if devid_temp in deviceid_replace:
            lines[my_id]['deviceid'] = deviceid_replace[devid_temp]

        # misc
        assert lines[my_id]['userid'] == "MITDELHI"

        for varname in ['country_code', 'full_address', 'hash_location', 'hash_precision', 'userid', 'timezone']:
            del lines[my_id][varname]

        # add UIDP from
        deviceid = lines[my_id]['deviceid']
        if deviceid in uidps:
            lines[my_id]['uidp'] = uidps[deviceid]
            lines[my_id]['name'] = names[deviceid]
        else:
            lines[my_id]['uidp'] = ''
            lines[my_id]['name'] = ''

        # accuracy minimum 10
        lines[my_id]['accuracy_level'] = max(10, int(float(lines[my_id]['accuracy_level'])))

        # event time
        temp = lines[my_id]['updated_time']
        update_time = datetime.strptime(temp[:-15], '%a %b %d %Y %H:%M:%S')
        if update_time.year != 2017:
            ids_to_delete.append(my_id)
        else:
            try:
                assert temp[-15:] in [' GMT+0000 (GMT)', ' GMT+0100 (BST)']
            except Exception as e:
                print(temp)
                print(my_id)
                raise e
            lines[my_id]['event_time'] = datetime.strptime(temp[:-15], time_format) \
                                         + timedelta(hours=5.5)
            # print(lines[my_id]['event_time'])
            del lines[my_id]['updated_time']

    # delete updates from other years
    for my_id in ids_to_delete:
        print("deleted 1 line because of year != 2017")
        del lines[my_id]

    print("finished processing lines")

    """
    sort by device id and time
    lines_sorted is a dictionary with deviceid and a list of points (sorted by time)
    """
    lines_list = [lines[my_id] for my_id in lines]
    deviceid_set = set([lines[my_id]['deviceid'] for my_id in lines])

    lines_by_devid = {}
    for line in lines_list:
        devid = line['deviceid']

        if devid not in lines_by_devid:
            lines_by_devid[devid] = [line]
        else:
            lines_by_devid[devid] += [line]

    print("finished sorting points")

    # return deviceid_set, lines_sorted
    return deviceid_set, lines_by_devid


def read_dump_deviceid(ROOT_PATH, line_buffer_size=500000):
    """
    Reads all data from the folder PATHIN, writes it to files by deviceid
    :param PATHIN:
    :param PATHOUT:
    :param singledate:
    :return:
    """

    print("BTR CODING - Reading ALL data")

    # paths
    pathsin = [ROOT_PATH + 'data/rawbtr/2017-02/',
               ROOT_PATH + 'data/rawbtr/2017-03/',
               ROOT_PATH + 'data/rawbtr/2017-04/',
               ROOT_PATH + 'data/rawbtr/2017-05/',
               ROOT_PATH + 'data/rawbtr/2017-06/',
               ROOT_PATH + 'data/rawbtr/2017-07/',
               ROOT_PATH + 'data/rawbtr/2017-08/',
               ROOT_PATH + 'data/rawbtr/2017-09/']

    PATHOUT = ROOT_PATH + "data/coded_deviceid/"

    # read existing uidp - deviceid correspondence
    temp = ut.csv2dict(ROOT_PATH + 'data/crosswalks/crosswalk_uidp_deviceid.csv', key_name='deviceid')
    names = {}
    uidps = {}
    for devid in temp:
        uidps[devid] = temp[devid]['uidp']
        names[devid] = temp[devid]['name']

    # read list of NEW deviceids to replace with OLD ones
    temp = ut.csv2dict(ROOT_PATH + 'data/crosswalks/deviceid_replacements.csv', key_name='deviceid_new')
    deviceid_replace = {}
    for deviceid_new in temp:
        deviceid_replace[deviceid_new] = temp[deviceid_new]['deviceid_old']

    file_list_to_read = []
    for PATHIN in pathsin:
        file_list_to_read += glob.glob(PATHIN + "*.csv")

    # chunks of data (multiple files)
    temp_data_chunk = {}
    reach_limit = False

    keys_to_delete = ["added_day", "added_hour", "added_min", "added_month", "added_sec", "added_time", "added_year",
                      "field_activity", "field_activityconf", "flag_archive", "from_time", "large_gap", "loc_altitude",
                      "loc_bearing", "loc_provider", "loc_speed", "loc_temp", "loc_wifi", "mode", "reverse_geo_status",
                      "time_gap_from_prev", "to_time", "version_name", "loc_humidity", "loc_state"]

    for idx, file in enumerate(file_list_to_read):
        file_time = file[-17:-4]

        """ ? """
        # last file
        if idx + 1 == len(file_list_to_read):
            reach_limit = True

        if not reach_limit:
            # clean rogue apostrophes
            clean_apostrophe_file(file, replace=True)

            # temp_data_chunk += ut.csv2dict(file)
            temp_data_chunk.update(ut.csv2dict(file, key_name='id'))

            print("BTR CODING - File " + file_time + " read. Length " + "{:,}".format(len(temp_data_chunk)))
            if len(temp_data_chunk) >= line_buffer_size:
                reach_limit = True

        else:
            try:
                # clean rogue apostrophes
                clean_apostrophe_file(file, replace=True)

                # temp_data_chunk += ut.csv2dict(file)
                temp_data_chunk.update(ut.csv2dict(file, key_name='id'))

                print("BTR CODING - File " + file_time + " read. Length " + "{:,}".format(len(temp_data_chunk)))
                print('BTR CODING - PROCESSING - !')

                # clean data
                deviceid_set, lines_sorted = code_btr_data(temp_data_chunk, uidps, names, deviceid_replace, debug=True)
                print('BTR CODING - finished coding data')

                # delete useless keys
                print("deleting useless keys - ")
                for devid in lines_sorted:
                    for line in lines_sorted[devid]:
                        for my_key in keys_to_delete:
                            if my_key in line:
                                del line[my_key]

                # append to files by deviceid
                for deviceid in deviceid_set:
                    print('BTR CODING - ' + file_time + " Processing deviceid = " + deviceid)

                    # current files paths
                    coded_btr_path = PATHOUT + 'coded_btr_' + deviceid + '.csv'

                    """ Save """
                    result = ut.dictlist2csv(lines_sorted[deviceid], coded_btr_path, append=True)
                    if result:
                        print('BTR CODING - ' + file_time + ' Success. Wrote full data to file: ' + coded_btr_path)
                    else:
                        print('BTR CODING - ' + file_time + ' Failed. Full data. Stopping now.')
                        sys.exit(-1)

            except Exception as e:
                print('BTR CODING - ' + file_time + " Error while loading BTR data file " + file)
                raise e

            # reset for next chunk
            temp_data_chunk = {}
            lines_sorted = None
            gc.collect()
            reach_limit = False

    return True


def clean_sort_deviceid(ROOT_PATH):

    # get list of files corresponding to deviceid
    file_list = glob.glob(ROOT_PATH + "data/coded_deviceid/*.csv")

    for file in sorted(file_list):
        devid = file[-13:-4]

        # read
        lines = ut.csv2dict(file)
        assert lines[0]['deviceid'] == devid
        print("BTR CODING - Processing deviceid " + devid)

        """ Sort by date """
        lines = sorted(lines, key=itemgetter('event_time'))

        """ Remove duplicate dates """
        event_time_set = set()
        temp = []
        deleted_lines = 0
        for line in lines:
            # try: del line['loc_humidity']
            # except: pass
            # try: del line['loc_state']
            # except: pass

            if line['event_time'] not in event_time_set:
                temp.append(line)
                event_time_set.add(line['event_time'])
            else:
                deleted_lines += 1

        lines = temp
        print('BTR CODING - ' + devid + " deleted " + str(deleted_lines) + " lines. (timestamp)")

        """ check id unique """
        idset = [line['id'] for line in lines]
        assert len(set(idset)) == len(idset)

        # only non-0 lat & long?
        # lines = [line for line in lines if line['latitude'] != '0']

        """ Save """
        result = ut.dictlist2csv(lines, file)
        if result:
            print('BTR CODING - ' + devid + ' Success. Wrote full data to file. ')
        else:
            print('BTR CODING - ' + devid + ' Failed. Full data. Stopping now.')
            sys.exit(-1)

    return True


def crosswalk_deviceid(ROOT_PATH):
    """
    create crosswalk with device ID and device name
    :param ROOT_PATH:
    :return:
    """
    import glob

    PATHIN = ROOT_PATH + "data/coded_deviceid/*.csv"
    PATHOUT = ROOT_PATH + "data/crosswalks/crosswalk_deviceid_devicename.csv"

    try:
        file_list_to_read = glob.glob(PATHIN)

        dictlist = []
        i = 0
        for i, file in enumerate(file_list_to_read):
            print(str(i) + "/" + str(len(file_list_to_read)))
            existing_data = ut.csv2dict(file)
            dictlist.append({
                'deviceid': existing_data[0]['deviceid'],
                'device_name': existing_data[0]['device_name'],
                'device_model': existing_data[0]['device_model']
            })

        # wite back to file
        result = ut.dictlist2csv(dictlist, PATHOUT)
        if result:
            print('Success writing crosswalk device ID - device name to file')
        else:
            print('Failed writing crosswalk device ID - device name to file')
            raise ArithmeticError

    except Exception as e:
        print("Error writing crosswalk device ID - device name to file")
        raise e


def crosswalk_check(ROOT_PATH):
    """
    Checks whether all the device IDs in the crosswalk_uidp_deviceid.csv EXIST in the data
    """
    import json
    from datetime import datetime

    now = str(datetime.now())[0:19].replace(" ", "_").replace(":", "")

    # refresh list of device ids
    # crosswalk_deviceid(ROOT_PATH)

    try:
        # read list of device ids
        lines1 = ut.csv2dict(ROOT_PATH + "data/raw_gps_other/crosswalk_deviceid_devicename.csv", key_name='deviceid')
        devidlist = list(lines1.keys())

        # read list of uidps and device ids
        lines2 = ut.csv2dict(ROOT_PATH + "data/coded_cto/crosswalk_uidp_deviceid.csv", key_name='deviceid')
        crosswalk = list(lines2.keys())

        # not so interesting
        print("In BTR data but not in survey")
        extra_devids = set(devidlist) - set(crosswalk)
        print(len(extra_devids))
        print(extra_devids)

        # important
        missing_devids = set(crosswalk) - set(devidlist)
        print("In survey but not in BTR data")
        print(len(missing_devids))
        print(missing_devids)

    except Exception as e:
        print("Error when checking deviceids in crosswalk")
        raise e

