print("importing packages \n")
import numpy as np
from scipy.spatial import cKDTree
import matplotlib as mpl
import pandas as pd
import geopandas as gpd
import os
from operator import itemgetter
from shapely.geometry import Polygon, shape, Point
import shapely.wkt
import time
import pathos.multiprocessing as mp
from pathos.parallel import stats
import sqlalchemy
from sqlalchemy import create_engine
import psycopg2
import csv
from scipy import spatial

#make postgres connection
host = 
dbase = 
user = 
pwd = 
con = create_engine('postgresql://{user}:{pwd}@{host}:5432/{dbase}'.format(host=host, dbase=dbase, user=user, pwd=pwd), echo=False)
connection = con.raw_connection()

# -----------------------------------------------------------------------------------------------------------

#                   DEFINE FUNCTIONS

# -----------------------------------------------------------------------------------------------------------
#get inputs function for debugging (debug=True), or to run for real (debug=False)
#get inputs function for debugging (debug=True), or to run for real (debug=False)
def get_inputs(debug=False):
    if debug == True:
        sql_bldgs = "select bid, the_geom_4326, the_point_4326 from la100.building_geoms limit 1000"#limit 1000
    else:
        sql_bldgs = "select bid, the_geom_4326, the_point_4326 from la100.building_geoms"
    poly_bldgs = gpd.GeoDataFrame.from_postgis(sql_bldgs, con, geom_col='the_point_4326')
    poly_bldgs = poly_bldgs.drop('the_geom_4326', axis=1)
    poly_bldgs = poly_bldgs.to_crs(epsg=3310)
    poly_bldgs = poly_bldgs.rename(columns={'the_point_4326': 'point_3310'}).set_geometry('point_3310')
    print("poly_bldgs: \n", poly_bldgs.head())

    if debug == True:
        sql_lann = "select bid, n_bid, dist_m from kwaechte.la_nn_all where bid in (select bid from la100.building_geoms limit 1000)"
    else:
        sql_lann = "select bid, n_bid, dist_m from kwaechte.la_nn_all"
    la_nn = pd.read_sql(sql_lann, con)
    print("la_nearest neighbors: \n", la_nn.head())

    if debug == True:
        sides = pd.read_csv('/projects/kwaechte/sides.csv')
        sides_test = sides.loc[sides.bid.isin(poly_bldgs.bid.values)].copy()
        crs = {'init': 'epsg:3310'}
        geometry = sides_test['point_3310'].map(shapely.wkt.loads)
        sides_test = sides_test.drop('point_3310', axis=1)
        sides_test = sides_test.drop('Unnamed: 0', axis=1)
        sides = gpd.GeoDataFrame(sides_test, crs=crs, geometry=geometry)
        sides['x'] = sides.geometry.x
        sides['y'] = sides.geometry.y
    else:
        sides = pd.read_csv('/projects/kwaechte/sides.csv')
        crs = {'init': 'epsg:3310'}
        geometry = sides['point_3310'].map(shapely.wkt.loads)
        sides = sides.drop('point_3310', axis=1)
        sides = sides.drop('Unnamed: 0', axis=1)
        sides = gpd.GeoDataFrame(sides, crs=crs, geometry=geometry)
        sides['x'] = sides.geometry.x
        sides['y'] = sides.geometry.y

    return la_nn, poly_bldgs, sides

#MULTIPROCESSING FUNCTIONS
def mp_main_worker(bids_split, cores, sides_split, bid_sides, la_nn_split):
    global BIDS_SPLIT
    BIDS_SPLIT = bids_split
    global SIDES_SPLIT
    SIDES_SPLIT = sides_split
    global BID_SIDES
    BID_SIDES = bid_sides
    global LA_NN_SPLIT
    LA_NN_SPLIT = la_nn_split
    print('init worker...')

def mp_process_worker(_id):
    # call in global vars
    sides = SIDES_SPLIT[_id].copy().reset_index()
    print("SIDES_SPLIT COUNT:", len(sides))
    bids = BIDS_SPLIT[_id]
    print("BIDS_SPLIT COUNT:", len(bids))
    bid_sides = BID_SIDES[_id].copy().reset_index()
    print("BID_SIDES COUNT:", len(bid_sides))
    la_nn = LA_NN_SPLIT[_id]
    print("LA_NN_SPLIT COUNT:", len(la_nn))

    # make interval list to print progress in terminal
    print_identifier = np.arange(0, len(sides), 1000)
    df = pd.DataFrame({'bid': [], 'uid': [], 'distance': [], 'neighbor': []})
    print("ENTERING FOR LOOP")

    #main work part:
    try:
        sides['processor_id'] = _id

        for bid in bids:

            try:
                # make df copy of building sides where self_bid != other_bids
                same_sides = bid_sides.loc[bid_sides.bid == bid].copy() ## pull from all of the sides where bid is not the same as loop's bid
                               nn_side_ids = la_nn.loc[la_nn.bid == bid].n_bid.values
                other_sides = sides.loc[sides.bid.isin(nn_side_ids)].copy().reset_index()

                for i, row in same_sides.iterrows():
                    nA = np.array(pd.DataFrame({'x': [row.x], 'y': [row.y]}))#, 'uid': [row.uid], 'bid': [row.bid]}))
                    nB = np.array(list(zip(other_sides.x, other_sides.y)))
                    dist, idx = spatial.KDTree(nB).query(nA)
#                     print(row.uid, ": dist:", dist[:], 'indx: ', other_sides.loc[idx[0], 'uid'], '; id:', _id)
                    df = df.append({'bid': row.bid, 'uid': row.uid, 'distance': dist[0], 'neighbor': other_sides.loc[idx[0]]}, ignore_index=True)
                    if i in print_identifier:
                        t_fcnstart = time.time()
                        checkpoint1 = time.ctime(int(t_fcnstart))
                        print(_id, ": ", checkpoint1, "// Place in df: ", i, "//", row.bid, dist)

            except Exception as ee:
                print(_id, bid, "|| ERROR(ee):", ee)

        return df

    except Exception as e:
        print("ERROR(e):", e)


def mp_func(sides, cores): #, bid_sides, la_nn_split, sides_split
    ids = np.arange(0, cores)
    # unique list of bids to chunk by
    bids = sides.bid.unique()
    print("Unsplit BIDs count:", len(bids))
    bids_split = np.array_split(bids, cores)

    #initialize lists for variable lists
    bid_sides = []
    la_nn_split = []
    sides_split = [] #list of lists for each building for which sides
    # make all the lists of lists
    for bids in bids_split:
        print("bids len:", len(bids))
        la_nn_split.append(la_nn.loc[la_nn.bid.isin(bids)].copy())
        nn_side_ids = la_nn.loc[la_nn.bid.isin(bids)].n_bid.values
        sides_split.append(sides.loc[sides.bid.isin(nn_side_ids)].copy())
        bid_sides.append(sides.loc[sides.bid.isin(bids)].copy())

    print("MOVING INTO MAIN PROCESS \n")
    main_worker_inputs = (bids_split, cores, sides_split, bid_sides, la_nn_split)
    pool = mp.Pool(cores, mp_main_worker, main_worker_inputs)
    results = pool.imap_unordered(mp_process_worker, ids)
    pool.close()
    pool.join()
    results_df = []
    for result in results:
        print("RESULTS LEN:", len(result), "\n", result.tail())
        results_df.append(result)
    offsets = pd.concat(results_df)
    print("\nRESULTS COMBINED\n")
    return offsets
    
# -----------------------------------------------------------------------------------------------------------

#                   CALL FUNCTIONS

# -----------------------------------------------------------------------------------------------------------
# get inputs
datastart = time.ctime(int(time.time()))
print("BRING IN DATA:", datastart)
la_nn, poly_bldgs, sides = get_inputs(debug=False)
dataend = time.ctime(int(time.time()))
print("BRING IN DATA COMPLETE:", dataend)
print('---------------------------------------------------------\n')

#start timing for mp_func
t0 = time.time()
start = time.ctime(int(t0))
print("\nStart time: ", str(start))

#### EXECUTE MULTIPROCESSING FUNCTION ###
coresnum = 28
offsets = mp_func(sides, coresnum) #, bid_sides, la_nn_split, sides_split

t1 = time.time()
end = time.ctime(int(t1))
print("\nEnd time: ", str(end))
print("\ntotal time elapsed for my_func: ", t1-t0, " seconds")

# -----------------------------------------------------------------------------------------------------------

#                   EXPORT RESULTS

# -----------------------------------------------------------------------------------------------------------
print('------------------EXPORT RESULTS--------------------------\n')
exportstart = time.ctime(int(time.time()))
print("EXPORTING DATA TO POSTGRES:", exportstart)
name = 'bldg_offsets'
schema = 'kwaechte'
offsets.to_sql(name=name, con=con, schema=schema, if_exists='replace', index=False, chunksize=20000) #dtype = {"bid":"numeric", "uid":"float", "distance":"float", "nearest_idx":"float"}
# con.close()
exportend = time.ctime(int(time.time()))
print("EXPORTED TO POSTGRES:", exportend)

# -----------------------------------------------------------------------------------------------------------

#                   CHECK RESULTS

# -----------------------------------------------------------------------------------------------------------
print('------------------CHECK RESULTS--------------------------\n')
offsets.distance.describe()
checking = offsets.groupby('bid').size().reset_index()
print("\nBUILDINGS WITH NOT 4 SIDES:", checking.loc[checking[0]!=4])
