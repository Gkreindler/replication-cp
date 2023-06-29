import csv
# import xlrd
from scipy.interpolate import interp1d
import numpy as np
from math import floor
from csv import DictReader
import os
import codecs
from pprint import pprint


def bool2int(my_dict):
    for key in my_dict.keys():
        if type(my_dict[key]) is bool:
            my_dict[key] = int(my_dict[key])
    return my_dict


def hf_hhmmss(tt):
    if len(tt) == 7:
        tt = tt.zfill(8)
    # compute fractional hour given 5 digit string hh:mm
    assert tt[2] == ":"
    hh = int(tt[0:2])
    mm = int(tt[3:5])
    return float(hh) + float(mm) / 60


def hf_hhmm(tt):
    if len(tt) == 4:
        tt = tt.zfill(5)
    # compute fractional hour given 5 digit string hh:mm
    assert tt[2] == ":"
    hh = int(tt[0:2])
    mm = int(tt[3:5])
    return float(hh) + float(mm) / 60


def get_colors():
    """
    Load Color Map
    """
    # colormapfile = "colormap_plasma.txt"
    colormapfile = "gps-maps/colormap_fire.txt"
    with open(colormapfile, 'r') as csvfile:
        colormap = list(csv.reader(csvfile, delimiter='\t'))

    # for fire
    colormap = [[int(l) / 255 for l in a] for a in colormap]

    # for plasma
    # colormap.reverse()

    # set up colormap by time of day
    xcm = np.linspace(0, 60 * 24 - 1, num=len(colormap), endpoint=True)
    cm_r = list(a[0] for a in colormap)
    cm_g = list(a[1] for a in colormap)
    cm_b = list(a[2] for a in colormap)

    f_r = interp1d(xcm, cm_r)
    f_g = interp1d(xcm, cm_g)
    f_b = interp1d(xcm, cm_b)

    hhmm = np.linspace(0, 60 * 24 - 1, num=60 * 24, endpoint=True)
    hhmm_r = f_r(hhmm)
    hhmm_g = f_g(hhmm)
    hhmm_b = f_b(hhmm)

    return hhmm_r, hhmm_g, hhmm_b


def clamp(x):
    return max(0, min(x, 255))


def range_color(x):
    """
    x is between 0 and 1
    :param x:
    :return:
    """
    hhmm_r, hhmm_g, hhmm_b = get_colors()
    try:
        # clamp
        x = min(1, max(0, x))
        x = floor(x * (24 * 60 - 1))
        r = floor(hhmm_r[x] * 255)
        g = floor(hhmm_g[x] * 255)
        b = floor(hhmm_b[x] * 255)

        return  "#{0:02x}{1:02x}{2:02x}".format(clamp(r), clamp(g), clamp(b))

    except Exception as e:
        print("failed prepping color for x =")
        print(x)
        raise e


def time_color(hour, minute, hhmm_r, hhmm_g, hhmm_b):
    try:
        # line color
        daymin = int(hour) * 60 + int(minute)
        # print(daymin)
        r = floor(hhmm_r[daymin] * 255)
        g = floor(hhmm_g[daymin] * 255)
        b = floor(hhmm_b[daymin] * 255)

        pointcolor = "#{0:02x}{1:02x}{2:02x}".format(clamp(r), clamp(g), clamp(b))

        return pointcolor
    except Exception as e:
        print(hour)
        print(minute)
        raise e


def xls2dict(file_name, sheet_index=0, key_name='', first_row=0):
    book = xlrd.open_workbook(filename=file_name)
    sheet = book.sheet_by_index(sheet_index)

    my_dict = [{
        sheet.cell_value(first_row, j): sheet.cell_value(i, j)
        for j in range(sheet.ncols)
    } for i in range(first_row + 1, sheet.nrows)]

    if key_name:
        try:
            assert len(set([line[key_name] for line in my_dict])) == len(my_dict)
        except Exception as e:
            print("key " + key_name + " not unique ")
            pprint(my_dict)
            raise e
        my_dict = {line[key_name]: line for line in my_dict}

    return my_dict


def csv2dict(pathin, key_name=''):
    """
    Self explanatory
    :param pathin:
    :param key_name:
    :return:
    """
    dicts = []

    if not os.path.isfile(pathin):
        return []

    with codecs.open(pathin, 'r', encoding='utf-8-sig') as csvfile:
        dicts = [d for d in DictReader(csvfile)]

    if key_name:
        try:
            mydict = {}
            for line in dicts:
                # make sure that key_name is indeed a unique key within the file
                assert line[key_name] not in mydict
                mydict[line[key_name]] = line
        except Exception as e:
            print('duplicate keys (' + key_name + ' in file ' + pathin)
            print(line)
            raise e

        return mydict
    else:
        return dicts


def csv2dict1horiz(pathin):
    """
    Reads a single dictionary written in a csv with keys on first column and values in second
    :param pathin:
    :return:
    """
    result = {}
    with codecs.open(pathin, 'r') as csvfile:
        lines = list(csv.reader(csvfile, delimiter=','))
        for line in lines:
            try:
                assert line[0] not in result
            except Exception as e:
                print("key " + str(line[0]) + " already read.")
                raise e

            result[line[0]] = line[1]
    return result


def csv_header():
    return ['id', 'idx', 'name', 'deviceid', 'uidp', 'battery', 'location_code',
            'date', 'time', 'time2prev', 'time2next',
            'latitude', 'longitude', 'dist2prev', 'dist2next', 'angle', 'accuracy_level',
            'trip', 'drop', 'drop_type',
            'score_before', 'score_after', 'markercolor', 'linecolor']


def csv2pointlist(pathin, silent=False):
    """
    Read all csv points files (taking the latest version of the manually edited csv)
    :param pathin: input folder
    :return:
    """
    import glob
    import re
    myheader = csv_header()

    """
    Decide which csv filesto read
    """
    file_list_v0 = glob.glob(pathin + "*.csv")
    file_dict = {}

    # convert list of files to dictionary by date, with a dictionary of files for each date (by modified date)
    for file in file_list_v0:
        p = re.compile("2017\-[0-9]{2}\-[0-9]{2}")
        file_date = p.findall(file)[0]

        stub = file[file.find("\\points_")+8:].replace(".csv", "")

        # add timestamp if necessary
        p = re.compile("_2017[0-9]{10}")
        has_timestamp = p.findall(stub)
        if not has_timestamp:
            modified_date = "20160101000000"
            stub += "_" + modified_date
        else:
            modified_date = has_timestamp[0][1:]

        if file_date not in file_dict:
            file_dict[file_date] = {}

        file_dict[file_date][modified_date] = {'stub': stub, 'file': file, 'date': file_date}

    # pick only the latest modified date for each date - still dictionary
    file_list_to_read = []
    for file_date in file_dict:
        mod_dates = list(file_dict[file_date].keys())
        last_mod_date = max(mod_dates)
        if not silent:
            print("Reading file " + file_dict[file_date][last_mod_date]['stub'] + " out of " + str(len(mod_dates)) + " files.")
        file_list_to_read.append(file_dict[file_date][last_mod_date])

    """
    read csv files
    dict[date] = (list of points)
    """
    if not silent:
        print(file_list_to_read)
    allpoints = {}
    for file in file_list_to_read:
        # allpoints[file['date']] = csv2dict(file['file'], key_name=u'id')
        allpoints[file['date']] = csv2dict(file['file'])

    return allpoints


def pointlist2csv(points_list, pathout, drop_color=True):
    myheader = csv_header()

    hhmm_r, hhmm_g, hhmm_b = get_colors()

    dictlist = []
    for point in points_list:
        date, time = str(point.timestamp).replace("+05:30", "").split(" ")
        hour = int(time[0:2])
        minute = int(time[3:5])

        mydict = {'id': point.pid,
                  'idx': point.idx,
                  'name': point.name,
                  'uidp': point.uidp,
                  'deviceid': point.deviceid,
                  'battery': point.battery_status,
                  'location_code': point.location_code,
                  'date': date,
                  'time': time,
                  'time2prev': point.time2prev,
                  'dist2prev': point.dist2prev,
                  'time2next': point.time2next,
                  'dist2next': point.dist2next,
                  'angle': point.angle,
                  'latitude': point.latitude,
                  'longitude': point.longitude,
                  'accuracy_level': point.accuracy_level,
                  'trip': int(point.trip),
                  'drop': int(point.drop),
                  'drop_only_color': int(point.drop_only_color),
                  'drop_type': point.drop_type,
                  'ignore': int(point.ignore),
                  'score_after': point.score_after,
                  'score_before': point.score_before,
                  'markercolor': time_color(hour, minute, hhmm_r, hhmm_g, hhmm_b),
                  'linecolor': time_color(hour, minute, hhmm_r, hhmm_g, hhmm_b)
                  }

        # do not actually drop, only color
        if point.drop_only_color and not drop_color:
            mydict['drop'] = 0
            mydict['markercolor'] = point.markercolor
            mydict['linecolor'] = point.markercolor

        dictlist.append(mydict)

    return dictlist2csv(dictlist, pathout, header=myheader)


def dictlist2csv(dictlist, pathout, header=None, append=False, na_replace='', cr_replace=False, debug=False):
    """
    Write a list of dictionaries to csv, (optionally) using an ordered list as headers.
    :param dictlist:
    :param pathout:
    :param header: when provided, it can be a subset of all keys in the dictionaries (union)
    :param append:
    :param na_replace:
    :param cr_replace: replace carriage return with space
    :param debug:
    :return:
    """
    from codecs import open

    # always write from scratch when file doesn't exist!
    if not os.path.isfile(pathout):
        append = False

    if append: file_mode = 'a'
    else:      file_mode = 'w'

    try:
        # If we're given a header, write in that order, and only those fields
        if header:
            rows_list = []
            if not append:
                rows_list.append(','.join(header) + '\n')
            for d in dictlist:
                row = [var2str4csv(d.get(field, na_replace)) for field in header]
                rows_list.append(','.join(row) + '\n')
            with open(pathout, file_mode, encoding='utf-8') as myfile:
                myfile.write(''.join(rows_list))

        else:  # If there is no header, write sorted keys (union)
            if debug:
                print("dictlist2csv: no header provided, using union of all headers, sorted alphabetically")
            keyset = set()
            for d in dictlist:
                keyset = set.union(keyset, set(d.keys()))
            keyset = sorted(list(keyset))

            if debug:
                print("dictlist2csv: finished preparing header")

            rows_list = []
            # write header
            if not append:
                rows_list.append(','.join([var2str4csv(k) for k in keyset]) + '\n')
            # write rows
            i = 0
            for d in dictlist:
                i += 1
                if i % 10000 == 1 and debug:
                    print("writing line + " + str(i))
                rows_list.append(','.join([var2str4csv(d.get(k, ''), cr_replace=cr_replace) for k in keyset]) + '\n')

            if debug:
                print("dictlist2csv: finished writing to variable, now writing to file:")

            with open(pathout, file_mode, encoding='utf-8') as myfile:
                myfile.write(''.join(rows_list))

        return True
    except Exception as e:
        print("dictlist2csv: Error printing to file " + pathout)
        print(e)
        return False


def json2dict(pathin, debug=False):
    from json import load
    from pprint import pprint
    from os.path import isfile

    if not isfile(pathin):
        return []

    with open(pathin) as data_file:
        data = load(data_file)

    if debug:
        pprint(data)

    return data


def stuff2json(dictlist, pathout):
    """
    write list of dictionaries to json file
    :param pathout:
    :return:
    """
    from json import dump
    if not dictlist:
        print("Empty json. No log file written.")
        return True
    try:
        with open(pathout, 'w') as fout:
            dump(dictlist, fout)
            print("Success saving json to file " + pathout)
            return True
    except Exception as e:
        print("error writing json to file " + pathout)
        print(e)
        return False


def var2str4csv(variable, cr_replace=False, time_format="%Y-%m-%d %H:%M:%S"):
    # "%d %b %Y %H:%M:%S"
    import datetime
    if variable == '':
        return ''
    if isinstance(variable, datetime.datetime):
        return "\"D " + variable.strftime(time_format) + "\""
    else:
        if cr_replace:
            return "\"" + str(variable).replace("\n", " ").replace("\r", " ") + "\""
        else:
            return "\"" + str(variable) + "\""


def str2hr(time_str, time_format):
    from time import strptime
    temp = strptime(time_str, time_format)
    hours = temp.tm_hour + temp.tm_min / 60 + temp.tm_sec / 3600
    return hours


def dictlist2dict(dictlist, field):
    """
    Take a list of dictionaries and make a dictionary of lists of dictionaries
    :param dictlist:
    :param field:
    :return:
    """


def str_color(r, g, b):
    return "#{0:02x}{1:02x}{2:02x}".format(clamp(r), clamp(g), clamp(b))


def clamp(x):
    return max(0, min(x, 255))


