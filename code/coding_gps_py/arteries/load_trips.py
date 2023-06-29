import json
import random
from pprint import pprint
from datetime import datetime
from geopy.distance import vincenty
import os
import utils as ut
import numpy as np
# import utils_geometry as ug
import googlemaps
import math
import matplotlib.pyplot as plt
import shapely.geometry as shgeo
import fiona
from pprint import pprint


def read_data_network():

    """ DEFINITIONS """
    precision = 5
    datapath = 'C:/bang_launch1/data/road_network/'  # Gabriel
    os.chdir(datapath)
    """ END """

    # read node locations
    with open('input/Nodes.geojson') as f:
        data_nodes = json.load(f)

    points = [{'lat': feature['geometry']['coordinates'][1],
               'long': feature['geometry']['coordinates'][0],
               'node_id': feature['properties']['name']} for feature in data_nodes['features']]

    # read undirected edges
    with open('input/Undirected-Edges.geojson') as g:
        data_edges = json.load(g)

    edges = []  # edge_id, origin_node_id, destination_node_id, curve_points
    errors = []

    for idx, edge in enumerate(data_edges['features']):
        print("processing edge " + str(idx))

        ori_lat = edge['geometry']['coordinates'][0][1]  # get latitude from first set of coordinates
        ori_long = edge['geometry']['coordinates'][0][0]  # get longitude from last set of coordinates
        dest_lat = edge['geometry']['coordinates'][-1][1]  # get latitude from last set of coordinates
        dest_long = edge['geometry']['coordinates'][-1][0]  # get longitude from last set of coordinates

        # distances to all points from originate point
        dist_meters_originate = [(point['node_id'], vincenty((ori_lat, ori_long), (point['lat'], point['long'])).meters)
                                 for point in points]
        # distances to all points from destination point
        dist_meters_destination = [(point['node_id'], vincenty((dest_lat, dest_long), (point['lat'], point['long'])).meters)
                                   for point in points]

        dist_meters_originate = sorted(dist_meters_originate, key=lambda x: x[1])  # sort ascending order to get closest point
        dist_meters_destination = sorted(dist_meters_destination, key=lambda x: x[1])
        closet_ori_point = dist_meters_originate[0]
        closet_dest_point = dist_meters_destination[0]

        # check uniqueness
        try:
            assert dist_meters_originate[1][1] > precision
        except Exception as e:
            print("Multiple close points Origin " + str(idx))
            pprint.pprint(edge)
            pprint.pprint(dist_meters_originate[0][0])
            pprint.pprint(dist_meters_originate[1][0])

        try:
            assert dist_meters_destination[1][1] > precision
        except Exception as e:
            print("Multiple close points Destination" + str(idx))
            pprint.pprint(edge)
            pprint.pprint(dist_meters_destination[0][0])
            pprint.pprint(dist_meters_destination[1][0])

        # add to list (assuming OK precision)
        if closet_ori_point[1] > precision or closet_dest_point[1] > precision:
            print("Error " + str(idx))
            pprint.pprint(edge)
            pprint.pprint(dist_meters_destination[0][0])
            errors.append({
                'name': edge['properties']['name'],
                'orig': closet_ori_point[0],
                'dest': closet_dest_point[0]
            })
        else:
            # invert lat-long
            temp = edge['geometry']['coordinates']
            temp = [(pt[1], pt[0]) for pt in temp]

            edges.append({
                'edge_route': str(temp),
                'name': edge['properties']['name'],
                'orig': closet_ori_point[0],
                'dest': closet_dest_point[0]
            })

    # write edges and points to output
    ut.dictlist2csv(points, 'network_clean/points.csv')
    ut.dictlist2csv(edges, 'network_clean/edges.csv')
    if errors:
        print("Errors found (edges that cannot be matched): " + str(len(errors)))
        ut.dictlist2csv(errors, 'network_clean/errors.csv')


def load_points_v0():
    """ DEFINITIONS """
    ROOT_PATH = 'C:/bang_launch1/'
    datapath = 'C:/bang_launch1/data/road_network/'
    os.chdir(ROOT_PATH)
    version = 'v4'
    chain_th = 15

    """ LOAD TRIPS """
    # load list of trips
    trip_sample_path = ROOT_PATH + 'data/road_network/input/trip_sample_full.csv'
    trip_sample = ut.csv2dict(trip_sample_path)

    trips_byuidp = {}

    # collect by uidp
    for trip in trip_sample:
        uidp = trip['uidp']
        date = trip['date']
        chai = trip['chain']
        if uidp in trips_byuidp:
            trips_byuidp[uidp].append(trip)
        else:
            trips_byuidp[uidp] = [trip]

    # by uidp: load all trips, keep only sample, add points
    allpoints_path = datapath + 'intermediate/all_points_v0.csv'
    try:
        os.remove(allpoints_path)
    except Exception as e:
        pass

    for idx, uidp in enumerate(sorted(trips_byuidp.keys())):
        allpoints = ''
        # list of date-chain
        date_chain_list = [(trip['date'], trip['chain']) for trip in trips_byuidp[uidp]]

        # load all chains for this uidp
        chains_data_file = ROOT_PATH + 'data/analysis/analysis_' + version + '_uidp/' + \
                           uidp + '/segs_chain_' + str(chain_th) + '.csv'
        chains_all = ut.csv2dict(chains_data_file)
        chains = [chain for chain in chains_all if (chain['date'], chain['chain']) in date_chain_list]
        # print("uidp: " + str(len(chains_all)) + ' chains, out of which ' + str(len(chains)) + ' in sample.')

        print(str(idx + 1) + ' of ' + str(len(trips_byuidp)))

        # load meta (segs)
        trip_data_file = ROOT_PATH + 'data/analysis/analysis_' + version + '_uidp/' + uidp + '/trips.csv'
        segs = ut.csv2dict(trip_data_file)
        # segs = [seg for seg in segs if seg['date'] == date and seg['type'] != 'end']

        """ Collect and add the points """
        for chain in chains:
            date = chain['date']
            # load points data - points_L0202093440_2017-02-18.csv
            points_data_file = ROOT_PATH + 'data/coded_uidp/' + \
                               uidp + '/points/points_' + \
                               uidp + '_' + date + '.csv'
            my_points = ut.csv2dict(points_data_file)
            my_points = [point for point in my_points]

            segs_list = eval(chain['segs'])
            my_segs = [seg for seg in segs if seg['date'] == date and seg['type'] != 'end']
            segs_dict = {int(seg['id']): seg for seg in my_segs}

            chain['points'] = []
            for idx in segs_list:
                seg = segs_dict[idx]
                if seg['type'] in ['jump', 'jumptrip', 'trip', 'gap']:
                    i_orig = int(seg['start'])
                    i_dest = int(segs_dict[idx + 1]['start'])
                    # chain['points'] += [(float(pt['latitude']), float(pt['longitude']), pt['time'])
                    #                     for pt in my_points if i_orig <= int(pt['idx']) <= i_dest
                    #                     and (pt['drop'] == '0' or i_orig == int(pt['idx']) or int(pt['idx']) == i_dest)]
                    temp = [pt['latitude'] + ',' + pt['longitude']
                            for pt in my_points if i_orig <= int(pt['idx']) <= i_dest
                            and (pt['drop'] == '0' or i_orig == int(pt['idx']) or int(pt['idx']) == i_dest)]
                    allpoints += '\n'.join(temp) + '\n'

        with open(allpoints_path, 'a') as myfile:
            myfile.write(allpoints)


def snap_trip(trip, gmaps_client=None):
    if not gmaps_client:
        mykey = ""
        gmaps_client = googlemaps.Client(key=mykey)

    # Split up the points into equal length lists of length <= 100
    rounds = []
    num_splits = math.ceil(len(trip) / 100)
    length = len(trip) // num_splits
    for i in range(num_splits - 1):
        rounds.append(trip[i * length: (i + 1) * length])
    rounds.append(trip[(num_splits - 1) * length:])

    result = []
    offset = 0
    for i, r in enumerate(rounds):
        # Run the algorithm
        print("number of points " + str(len(r)))
        res = gmaps_client.snap_to_roads(r, interpolate=True)
        # pprint(res)

        # Offset the originalIndices to compensate for splitting into rounds
        for x in res:
            if 'originalIndex' in x:
                x["originalIndex"] = x["originalIndex"] + offset
        offset += len(r)
        result.extend(res)

    print("Done. Snapped " + str(len(result)) + " points out of " + str(len(trip)) + " points.")
    return result


def snap_all():
    """ load trip points, snap with google, and  """
    ROOT_PATH = 'C:/bang_launch1/'
    datapath = 'C:/bang_launch1/data/road_network/'
    os.chdir(ROOT_PATH)
    version = 'v4'
    chain_th = 15

    random.seed(a=23492743)

    """ Google Maps API client """
    mykey = ""
    gmaps_client = googlemaps.Client(key=mykey)

    """ LOAD TRIPS """
    # load list of trips
    trip_sample_path = ROOT_PATH + 'data/road_network/input/trip_sample_full.csv'
    trip_sample = ut.csv2dict(trip_sample_path)

    trips_byuidp = {}

    # collect by uidp
    for trip in trip_sample:
        uidp = trip['uidp']
        date = trip['date']
        chai = trip['chain']
        if uidp in trips_byuidp:
            trips_byuidp[uidp].append(trip)
        else:
            trips_byuidp[uidp] = [trip]

    # for idx, uidp in enumerate(sorted(trips_byuidp.keys())):

    keys_list = sorted(trips_byuidp.keys())
    random.shuffle(keys_list)
    keys_list = keys_list[:100]
    all_trips = []
    all_snapd = []
    all_snapd_file = datapath + 'intermediate/all_snapd_temp.csv'
    all_trips_file = datapath + 'intermediate/all_trips_temp.csv'

    for idx, uidp in enumerate(keys_list):
        # list of date-chain
        date_chain_list = [(trip['date'], trip['chain']) for trip in trips_byuidp[uidp]]

        # todo: one per uidp for now
        date_chain_list = [random.choice(date_chain_list)]

        # load all chains for this uidp
        chains_data_file = ROOT_PATH + 'data/analysis/analysis_' + version + '_uidp/' + \
                           uidp + '/segs_chain_' + str(chain_th) + '.csv'
        chains_all = ut.csv2dict(chains_data_file)
        chains = [chain for chain in chains_all if (chain['date'], chain['chain']) in date_chain_list]
        # print("uidp: " + str(len(chains_all)) + ' chains, out of which ' + str(len(chains)) + ' in sample.')

        print(str(idx + 1) + ' of ' + str(len(trips_byuidp)))

        # load meta (segs)
        trip_data_file = ROOT_PATH + 'data/analysis/analysis_' + version + '_uidp/' + uidp + '/trips.csv'
        segs = ut.csv2dict(trip_data_file)
        # segs = [seg for seg in segs if seg['date'] == date and seg['type'] != 'end']

        """ Collect and add the points """
        for chain in chains:
            date = chain['date']
            # load points data - points_L0202093440_2017-02-18.csv
            points_data_file = ROOT_PATH + 'data/coded_uidp/' + \
                               uidp + '/points/points_' + \
                               uidp + '_' + date + '.csv'
            my_points = ut.csv2dict(points_data_file)
            my_points = [point for point in my_points]

            segs_list = eval(chain['segs'])
            my_segs = [seg for seg in segs if seg['date'] == date and seg['type'] != 'end']
            segs_dict = {int(seg['id']): seg for seg in my_segs}

            chain['points'] = []
            for idx in segs_list:
                seg = segs_dict[idx]
                if seg['type'] in ['jump', 'jumptrip', 'trip', 'gap']:
                    i_orig = int(seg['start'])
                    i_dest = int(segs_dict[idx + 1]['start'])
                    chain['points'] += [(float(pt['latitude']), float(pt['longitude']), pt['time'])
                                        for pt in my_points if i_orig <= int(pt['idx']) <= i_dest
                                        and (pt['drop'] == '0' or i_orig == int(pt['idx']) or int(pt['idx']) == i_dest)]

                    # temp = [pt['latitude'] + ',' + pt['longitude']
                    #         for pt in my_points if i_orig <= int(pt['idx']) <= i_dest
                    #         and (pt['drop'] == '0' or i_orig == int(pt['idx']) or int(pt['idx']) == i_dest)]
                    # allpoints += '\n'.join(temp) + '\n'

            """ Snap to Roads """
            trip_wo_time = [(x[0], x[1]) for x in chain['points']]
            trip_snap = snap_trip(trip_wo_time, gmaps_client)

            """ add """
            all_trips += [{
                'uidp': uidp,
                'date': date,
                'id': chain['chain'],
                'lat': x[0], 'lon': x[1]
                          } for x in chain['points']]
            all_snapd += [{
                              'uidp': uidp,
                              'date': date,
                              'chain': chain['chain'],
                              'lat': x['location']['latitude'], 'lon': x['location']['longitude'],
                              'oid': x.get('originalIndex', '')
                          } for x in trip_snap]

            """ plot """
            # trip_snap_x = [x['location']['latitude'] for x in trip_snap]
            # trip_snap_y = [x['location']['longitude'] for x in trip_snap]
            #
            # trip_x = [x[0] for x in chain['points']]
            # trip_y = [x[1] for x in chain['points']]
            #
            # fig, ax1 = plt.subplots()
            # ax1.plot(trip_snap_x, trip_snap_y, 'r', label='snap')  # , linewidth=1.0
            # ax1.plot(trip_x, trip_y, 'b', label='original')  # , linewidth=1.0
            # plt.show()

            # # Make google map to check quality
            # ug.heatmap_trip(overlap_trips=[list(zip(trip_snap_x, trip_snap_y))],
            #                 non_overlap_trips=[trip_wo_time], edge_alpha=1, edge_width=2, edge_color='red',
            #                 markers=True, file_out=datapath + 'temp/map_test.html')

    """ Save to file """
    ut.dictlist2csv(all_snapd, all_snapd_file)
    ut.dictlist2csv(all_trips, all_trips_file)


def full_trip_list():
    """ load trip points, snap with google, and  """
    ROOT_PATH = 'C:/bang_launch1/'
    datapath = 'C:/bang_launch1/data/road_network/'
    os.chdir(ROOT_PATH)
    version = 'v4'
    chain_th = 15

    # random.seed(a=23492743)

    # shape = fiona.open("my_shapefile.shp")
    shapes = fiona.open("C:/bang_launch1/data/road_network/qgis/road_segment_btm.shp")

    # check ids are unique
    ids = [shape['id'] for shape in shapes]
    assert(len(ids) == len(set(ids)))

    pprint(ids)

    for myshape in shapes:
        # print(shape.schema)
        # first feature of the shapefile
        # myshape = shape.next()  # 1: edges 50, 51
        # myshape = shape.next()  # 2: edges 100, 100
        # myshape = shape.next()  # 3: edges 106, 107
        pprint(myshape)  # (GeoJSON format)

        mypoly = shgeo.shape(myshape['geometry'])

        # point = shgeo.Point(77.585999000000001, 12.931659000000000)
        # # print(mypoly.contains(point))
        # print(mypoly.distance(point) * 111.325)
        # dsfdfsdsf

        """ LOAD TRIPS """
        # load list of trips
        trip_sample_path = ROOT_PATH + 'data/road_network/input/trip_sample_full.csv'
        trip_sample = ut.csv2dict(trip_sample_path)

        trips_byuidp = {}

        # collect by uidp
        for trip in trip_sample:
            uidp = trip['uidp']
            date = trip['date']
            chai = trip['chain']
            if uidp in trips_byuidp:
                trips_byuidp[uidp].append(trip)
            else:
                trips_byuidp[uidp] = [trip]

        # for idx, uidp in enumerate(sorted(trips_byuidp.keys())):

        keys_list = sorted(trips_byuidp.keys())
        random.shuffle(keys_list)
        # keys_list = keys_list[:100]
        # all_trips = []
        all_trips_btm = []
        # all_trips_file = datapath + 'intermediate/all_trips_btm.csv'
        all_trips_btm_file = datapath + 'btm_volumes/all_trips_btm_pts_' + myshape['id'] + '.csv'

        for idx, uidp in enumerate(keys_list):
            # list of date-chain
            date_chain_list = [(trip['date'], trip['chain']) for trip in trips_byuidp[uidp]]

            # todo: one per uidp for now
            # date_chain_list = [random.choice(date_chain_list)]

            # load all chains for this uidp
            chains_data_file = ROOT_PATH + 'data/analysis/analysis_' + version + '_uidp/' + \
                               uidp + '/segs_chain_' + str(chain_th) + '.csv'
            chains_all = ut.csv2dict(chains_data_file)
            chains = [chain for chain in chains_all if (chain['date'], chain['chain']) in date_chain_list]
            # print("uidp: " + str(len(chains_all)) + ' chains, out of which ' + str(len(chains)) + ' in sample.')

            print(str(idx + 1) + ' of ' + str(len(trips_byuidp)))

            # load meta (segs)
            trip_data_file = ROOT_PATH + 'data/analysis/analysis_' + version + '_uidp/' + uidp + '/trips.csv'
            segs = ut.csv2dict(trip_data_file)
            # segs = [seg for seg in segs if seg['date'] == date and seg['type'] != 'end']

            """ Collect and add the points """
            for chain in chains:
                print('.', end="", flush=True)
                date = chain['date']
                # load points data - points_L0202093440_2017-02-18.csv
                points_data_file = ROOT_PATH + 'data/coded_uidp/' + \
                                   uidp + '/points/points_' + \
                                   uidp + '_' + date + '.csv'
                my_points = ut.csv2dict(points_data_file)
                my_points = [point for point in my_points]

                segs_list = eval(chain['segs'])
                my_segs = [seg for seg in segs if seg['date'] == date and seg['type'] != 'end']
                segs_dict = {int(seg['id']): seg for seg in my_segs}

                chain['points'] = []
                for idx in segs_list:
                    seg = segs_dict[idx]
                    if seg['type'] in ['jump', 'jumptrip', 'trip', 'gap']:
                        i_orig = int(seg['start'])
                        i_dest = int(segs_dict[idx + 1]['start'])
                        chain['points'] += [{'lat': float(pt['latitude']), 'lng': float(pt['longitude']), 't': pt['time']}
                                            for pt in my_points if i_orig <= int(pt['idx']) <= i_dest
                                            and (pt['drop'] == '0' or i_orig == int(pt['idx']) or int(pt['idx']) == i_dest)]

                """ Check if intersects Polygon """
                n_inside = 0
                for point in chain['points']:
                    shpt = shgeo.Point(point['lng'], point['lat'])
                    point['c'] = int(mypoly.distance(shpt) * 111.325 < 0.1)
                    n_inside += point['c']

                # frac_btm = np.mean([point['c'] for point in chain['points']])
                # all_trips.append({
                #     'frac_btm': frac_btm,
                #     'date': date,
                #     'start': chain['points'][0]['t']
                # })

                """ PLOT """
                # trip_coords = [(point['lat'], point['lng']) for point in chain['points']]
                # trip_btm_co = [(point['lat'], point['lng']) for point in chain['points'] if point['c'] == 1]

                # Make google map to check quality
                if n_inside > 1:
                    # ug.heatmap_trip(overlap_trips=[trip_btm_co],
                    #                 non_overlap_trips=[trip_coords], edge_alpha=0.5, edge_width=5, edge_color='red',
                    #                 markers=True, file_out=datapath + 'temp/map_test.html')
                    # sfdsfds
                    trip_btm = [point for point in chain['points'] if point['c'] == 1]

                    trip_length = np.sum([vincenty((trip_btm[i  ]['lat'], trip_btm[i  ]['lng']),
                                                   (trip_btm[i+1]['lat'], trip_btm[i+1]['lng'])).meters
                                          for i in range(len(trip_btm) - 1)])

                    t0 = datetime.strptime(trip_btm[0]['t'], '%H:%M:%S')
                    t1 = datetime.strptime(trip_btm[-1]['t'], '%H:%M:%S')

                    all_trips_btm.append({
                        'date': date,
                        'start': trip_btm[0]['t'],
                        'dist': vincenty((trip_btm[0]['lat'], trip_btm[0]['lng']), (trip_btm[-1]['lat'], trip_btm[-1]['lng'])).kilometers,
                        'dur': (t1-t0).total_seconds() / 3600,
                        'dir_lat': (trip_btm[-1]['lat'] - trip_btm[0]['lat']) * 111.325,
                        'dir_lng': (trip_btm[-1]['lng'] - trip_btm[0]['lng']) * 111.325,
                        'length': trip_length
                    })
            print('.')
        """ Save to file """
        # ut.dictlist2csv(all_trips, all_trips_file)
        ut.dictlist2csv(all_trips_btm, all_trips_btm_file)


if __name__ == '__main__':
    read_data_network()
    snap_all()
    full_trip_list()
