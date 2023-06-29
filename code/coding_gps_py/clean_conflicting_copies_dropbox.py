import glob
import os
import utils as ut
from pprint import pprint


def clean_subfolder(folderpath, readonly=True):
    filelist = glob.glob(folderpath + '*\'s conflicted copy *')
    for file in filelist:
        print(file)
        if not readonly:
            os.remove(file)
            print("deleted")


def clean_tree(directory, readonly=True):
    dirlist = [x[0] for x in os.walk(directory)][0:]
    for dir_path in dirlist:
        clean_subfolder(dir_path + '/', readonly=readonly)


if __name__ == '__main__':
    readonly = False

    ROOT_PATH = 'C:/bang_launch1/'
    clean_tree(ROOT_PATH + 'data/', readonly=readonly)
    clean_tree(ROOT_PATH + 'data/coded_uidp/', readonly=readonly)
    clean_tree(ROOT_PATH + 'data/analysis/analysis_v3_uidp/', readonly=readonly)
    clean_tree(ROOT_PATH + 'do_not_touch_this_folder/', readonly=readonly)
    clean_tree(ROOT_PATH + 'do_not_touch_this_folder/treat_data/', readonly=readonly)

    # only works on GK's computer
    try:
        clean_tree('C:/Users/Gabriel K/Dropbox (MIT)/0Re Projects/bangaloretraffic/code/bang-treat/', readonly=readonly)
        clean_tree('C:/Users/Gabriel K/Dropbox/0Re Projects/b_traffic_pilot2/code/circle-detour/', readonly=readonly)
        clean_tree("C:/Users/Gabriel K/Dropbox/0Re Projects/bangaloretraffic/code/gps-trips-bang/", readonly=readonly)
        clean_tree("C:/Users/Gabriel K/Dropbox (MIT)/0Re Projects/bangaloretraffic/code/road_network/", readonly=readonly)

    except:
        print("error")





