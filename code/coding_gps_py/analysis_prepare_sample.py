import utils as ut
from pprint import pprint
from operator import itemgetter
import numpy as np
from geopy.distance import geodesic
import glob
import os


def rint(mystr):
    """ convert string to int after rounding """
    return int(np.round(float(mystr)))


# def check_quality(uidp, date):
#     qual = ut.csv2dict('C:/bang_launch1/data/monitoring/quality_gaps/quality_all_new.csv')
#     return [line['data_quality'] for line in qual if line['uidp'] == uidp and line['date'] == date]


def get_chain_info_sid(field, date, sids, structs, structs_address):
    values = []
    for sid in sids:
        key = date + '_' + str(sid)
        idx = structs_address[key]

        # check:
        assert structs[idx]['date'] == date
        assert structs[idx]['sid'] == str(sid)

        values.append(structs[idx][field])
    return values


def prepare_sample(ROOT_PATH, uidps, PATHIN_STUB, PATHOUT=None, chain_th=None, sample=''):
    """
    Load chains and other files and trip sample selection:
        -- quality in dt precision
        -- quality in terms of jumps
        -- outstaion
    Export working sample
    """
    assert chain_th in [15, 30, 60]

    # read confirmed coded hw locations
    confirmed_hw_file = ROOT_PATH + 'data/raw_gps_homework/hw_confirmed.csv'
    confirmed_hw = ut.csv2dict(confirmed_hw_file, key_name='uidp')

    all_chains = []
    for idx, uidp in enumerate(uidps):
        print(str(idx) + '/' + str(len(uidps)) + ". processing " + uidp)

        # load chains
        chain_file = PATHIN_STUB + uidp + '/segs_chain_' + str(chain_th) + '.csv'
        chains = ut.csv2dict(chain_file)
        try:
            assert chains
        except Exception as e:
            print("NO chain " + str(chain_th) + " for uidp " + uidp)

        # load structures
        struct_file = PATHIN_STUB + uidp + '/segs_struct.csv'
        structs = ut.csv2dict(struct_file)
        # address based on date and SID
        structs_address = {line['date'] + '_' + line['sid']: idx for idx, line in enumerate(structs)}
        assert len(structs_address) == len(structs)

        if structs:
            # check that chains IDs are included
            chain_id_key = 'chain_' + str(chain_th)
            assert chain_id_key in structs[0]

            """ load HW info """
            hw_locs = confirmed_hw.get(uidp, None)
            work2_pt = work3_pt = None
            if hw_locs:
                home_pt = eval(hw_locs['orig'])
                work_pt = eval(hw_locs['dest'])
                if hw_locs['dest2']:
                    work2_pt = eval(hw_locs['dest2'])
                if hw_locs['dest3']:
                    work3_pt = eval(hw_locs['dest3'])
            else:
                home_pt = work_pt = (0.0, 0.0)
            # locations = home_pt, work_pt, work2_pt, work3_pt

            # define tolerances hw, hw2, ww2
            if hw_locs:
                dist_hw = geodesic(home_pt, work_pt).meters
                dist_hw_tol = 250 + min(750, dist_hw / 5000 * 750)
                if work2_pt:
                    dist_hw2 = geodesic(home_pt, work2_pt).meters
                    dist_hw2_tol = 250 + min(750, dist_hw2 / 5000 * 750)

                    dist_ww2 = geodesic(work_pt, work2_pt).meters
                    dist_ww2_tol = 250 + min(750, dist_ww2 / 5000 * 750)
                else:
                    dist_hw2_tol = -1
                    dist_ww2_tol = -1
            else:
                dist_hw_tol = dist_hw2_tol = dist_ww2_tol = -1
            assert not work3_pt

            # load OD info todo: DO THIS WHEN READY

            # collect chains
            uidp_chains = []
            for idx_chain, chain in enumerate(chains):
                assert int(chain['trip']) in [0, 1]
                if int(chain['trip']) == 1:
                    """ seg ids """
                    date = chain['date']
                    sids = eval(chain['sids'])
                    chain_id = int(chain['chain'])

                    """ outstation """
                    if float(chain['pl']) > 0:
                        chain['inb_frac'] = float(chain['plb']) / float(chain['pl'])
                    else:
                        chain['inb_frac'] = 0

                    """ total duration at ALL locations """
                    sigs_durs = get_chain_info_sid('dur_bg_', date, sids, structs, structs_address)
                    sigs_type = get_chain_info_sid('trip', date, sids, structs, structs_address)
                    # take only locations
                    locs_durs = [float(item[0]) for item in zip(sigs_durs, sigs_type) if item[1] == '0']
                    chain['locs_dur_bg'] = sum(locs_durs) / 60.0  # in minutes for comparison with dur_bg_mm

                    """ quality: departure time precision """
                    chain['t_start_prec'] = get_chain_info_sid('t_start_prec', date, sids[0:1], structs, structs_address)
                    chain['t_start_prec'] = int(np.round(int(chain['t_start_prec'][0])/60))
                    chain['t_stop_prec'] = get_chain_info_sid('t_stop_prec', date, sids[0:1], structs, structs_address)
                    chain['t_stop_prec'] = int(np.round(int(chain['t_stop_prec'][0])/60))

                    sigs_start_prec = get_chain_info_sid('t_start_prec', date, sids, structs, structs_address)
                    sigs_stop_prec = get_chain_info_sid('t_stop_prec', date, sids, structs, structs_address)

                    # only for trips
                    trips_start_prec = [int(np.round(int(item[0]) / 60.0)) for item in zip(sigs_start_prec, sigs_type) if item[1] == '1']
                    trips_stop_prec = [int(np.round(int(item[0]) / 60.0)) for item in zip(sigs_stop_prec, sigs_type) if item[1] == '1']
                    chain['trips_start_prec'] = str(trips_start_prec)
                    chain['trips_stop_prec'] = str(trips_stop_prec)

                    chain['t_start_prec_max'] = max(trips_start_prec)
                    chain['t_stop_prec_max'] = max(trips_stop_prec)

                    """ quality: jumps """
                    # todo: idea: reasonable speed jumps are OK
                    # todo: get vector of jump length and duration
                    # todo: need to go to structs get vector of jump length and duration
                    # field = 'jump'
                    trips_dur_jump = get_chain_info_sid('dur_jump', date, sids, structs, structs_address)
                    trips_dur_jump = sum([float(item[0]) for item in zip(trips_dur_jump, sigs_type) if item[1] == '1'])
                    chain['trips_dur_jump'] = rint(trips_dur_jump)
                    trips_plb_jump = get_chain_info_sid('plb_jump', date, sids, structs, structs_address)
                    trips_plb_jump = sum([float(item[0]) for item in zip(trips_plb_jump, sigs_type) if item[1] == '1'])
                    chain['trips_plb_jump'] = rint(trips_plb_jump)

                    """ classifty chain endpoints with HW information """
                    if chain['center_orig_lat'] == '' or chain['center_orig_lon'] == '':
                        assert chain['chain'] == '0'
                        orig = (0.0, 0.0)
                    else:
                        orig = (float(chain['center_orig_lat']), float(chain['center_orig_lon']))

                    if chain['center_dest_lat'] == '' or chain['center_dest_lon'] == '':
                        assert idx_chain == len(chain) - 1
                        dest = (0.0, 0.0)
                    else:
                        dest = (float(chain['center_dest_lat']), float(chain['center_dest_lon']))

                    # tolerances
                    chain['dist_hw_tol'] = dist_hw_tol
                    chain['dist_hw2_tol'] = dist_hw2_tol
                    chain['dist_ww2_tol'] = dist_ww2_tol

                    # distances HW
                    chain['d_oh'] = geodesic(orig, home_pt).meters
                    chain['d_ow'] = geodesic(orig, work_pt).meters
                    chain['d_dh'] = geodesic(dest, home_pt).meters
                    chain['d_dw'] = geodesic(dest, work_pt).meters

                    # close (HW)
                    chain['oh'] = oh = geodesic(orig, home_pt).meters < dist_hw_tol
                    chain['ow'] = ow = geodesic(orig, work_pt).meters < dist_hw_tol
                    chain['dh'] = dh = geodesic(dest, home_pt).meters < dist_hw_tol
                    chain['dw'] = dw = geodesic(dest, work_pt).meters < dist_hw_tol

                    # chain['h'] = h = oh or dh
                    # chain['w'] = w = ow or dw
                    # chain['hw'] = hw = oh and dw
                    # chain['wh'] = wh = ow and dh

                    # distances and close (W2)
                    if work2_pt:
                        chain['d_ow2'] = geodesic(orig, work2_pt).meters
                        chain['d_dw2'] = geodesic(dest, work2_pt).meters

                        chain['ow2'] = ow2 = geodesic(orig, work2_pt).meters < dist_hw2_tol
                        chain['dw2'] = dw2 = geodesic(dest, work2_pt).meters < dist_hw2_tol

                    for field in ['oh', 'ow', 'dh', 'dw', 'ow2', 'dw2']:
                        if field in chain:
                            chain[field] = int(chain[field])

                    """ get OD information """
                    # todo: loop through chains and classify

                    """ done! add trips """
                    uidp_chains.append(chain)

            """ done! add trips """
            all_chains += uidp_chains

            """ write to uidp folder if not collecting ALL """
            if not PATHOUT:
                sample_file = PATHIN_STUB + uidp + '/sample_chains_' + str(chain_th) + '.csv'
                if ut.dictlist2csv(uidp_chains, sample_file): print("Success!")
                else: print("Failed writing to file")

    """ write ALL trip chains to PATHOUT """
    if PATHOUT:
        sample_file = PATHOUT + '/sample_' + sample + 'chains_' + str(chain_th) + '.csv'
        if ut.dictlist2csv(all_chains, sample_file): print("Success!")
        else: print("Failed writing to file")


def update_chain_in_struct(uidps, PATHIN_STUB, chain_th):
    """
    Update chain info in the struct files
    :chain_th: chain threshold (15, 30, 60)
    """

    assert chain_th in [15, 30, 60]

    for idx, uidp in enumerate(uidps):
        print(str(idx) + '/' + str(len(uidps)) + ". processing " + uidp)
        chain_file = PATHIN_STUB + uidp + '/segs_chain_' + str(chain_th) + '.csv'
        chains = ut.csv2dict(chain_file)

        struct_file = PATHIN_STUB + uidp + '/segs_struct.csv'
        structs = ut.csv2dict(struct_file)
        # address based on date and SID
        structs_address = {line['date'] + '_' + line['sid']: idx for idx, line in enumerate(structs)}
        assert len(structs_address) == len(structs)

        for chain in chains:
            date = chain['date']
            sids = eval(chain['sids'])
            chain_id = int(chain['chain'])

            for sid in sids:
                key = date + '_' + str(sid)
                idx = structs_address[key]

                # check:
                assert structs[idx]['date'] == date
                assert structs[idx]['sid'] == str(sid)

                # add chain ID information to struct
                structs[idx]['chain_' + str(chain_th)] = chain_id

        # write back to file
        struct_file = PATHIN_STUB + uidp + '/segs_struct.csv'
        ut.dictlist2csv(structs, struct_file)
    return True


def prepare_sample_complete(ROOT_PATH, chain_th):

    PATHIN_STUB = ROOT_PATH + 'data/coded_gps_byid/'
    PATHOUT = ROOT_PATH + 'data/coded_gps/'

    """ LOAD ALL UIDPS """
    ALL_UIDPS_PATH = ROOT_PATH + 'data/coded_cto/crosswalk_uidp_deviceid.csv'
    uidp_list = ut.csv2dict(ALL_UIDPS_PATH)
    uidps = [line['uidp'] for line in uidp_list if
                 os.path.isfile(ROOT_PATH + 'data/raw_gps/coded_btr_' + line['deviceid'] + '.csv')]
    print(len(uidps))

    prepare_sample(ROOT_PATH=ROOT_PATH, uidps=uidps, PATHIN_STUB=PATHIN_STUB, PATHOUT=PATHOUT,
                   chain_th=chain_th, sample='complete_')

    return True


def prepare_sample_experiment(ROOT_PATH, chain_th):
    PATHIN_STUB = ROOT_PATH + 'data/coded_gps_byid/'
    PATHOUT = ROOT_PATH + 'data/coded_gps/'

    # EXPERIMENT uidps
    TREAT_PATH = ROOT_PATH + 'data/treatment/treatment roster noPII.csv'
    treat_roster = ut.csv2dict(TREAT_PATH)
    uidps = [line['uidp'] for line in treat_roster if line['meeting'] == 'done']
    # uidps = ['L0202100848']

    prepare_sample(ROOT_PATH=ROOT_PATH, uidps=uidps, PATHIN_STUB=PATHIN_STUB, PATHOUT=PATHOUT,
                   chain_th=chain_th, sample='exp_')

    return True


def count_by_id(PATH_BASE):
    folder_list = [name for name in os.listdir(PATH_BASE) if os.path.isdir(os.path.join(PATH_BASE, name))]

    # folder_list = glob.glob(PATHOUT)

    n_files_list = []
    for uidp in folder_list:
        print(uidp)

        my_path = PATH_BASE + uidp + '/points/*.csv'
        n_files = len(glob.glob(my_path))

        n_files_list.append({
            'uidp': uidp,
            'n_files': n_files
        })

    return n_files_list


# if __name__ == '__main__':
#
#     ROOT_PATH = 'C:/bang_cp_paper/'
#     PATHIN_STUB = ROOT_PATH + 'data/coded_uidp/'
#     PATHOUT = ROOT_PATH + 'data/coded_gps/'
#
#     """ check """
#     # print(check_quality('L0207085326', '2017-02-14'))
#
#     """ count files """
#     # PATH_BASE = ROOT_PATH + 'data/coded_uidp/'
#     # n_files_list = count_by_id(PATH_BASE)
#     # ut.dictlist2csv(n_files_list, ROOT_PATH + 'data/analysis/analysis_' + version + '_other/n_files_list.csv')


