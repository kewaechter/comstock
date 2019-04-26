import numpy as np
import matplotlib as mpl
import pandas as pd
from scipy.spatial import cKDTree
import geopandas as gpd
import psycopg2
import os
from operator import itemgetter
from shapely.geometry import Polygon, shape, Point
import shapely.wkt
import time

#make pg connection
host = ''
dbase = ''
user = ''
pwd = ''
con = psycopg2.connect(database=dbase, user=user, password=pwd, host=host)
sql_bldgs = "select bid, the_geom_4326, the_point_4326 from la100.building_geoms limit 1000" #************

#load data as geodataframes
poly_bldgs = gpd.GeoDataFrame.from_postgis(sql_bldgs, con, geom_col='the_geom_4326')
poly_bldgs = poly_bldgs.drop('the_point_4326', axis=1)
print("Building polys CRS: ", poly_bldgs.crs, "\n")

t0 = time.time()
# reproject data into projected CRS
poly_bldgs = poly_bldgs.to_crs(epsg=3310)
print("New CRS: ",poly_bldgs.crs, "\n")

print("starting directional coordinates")
# make new dataframe with directional extreme coordinates
## blank df to population with loop
df = pd.DataFrame({'index': [], 'bid': [], 'north': [], 'east': [], 'south': [], 'west': []})

for i, row in poly_bldgs.iterrows():
    geom = row.the_geom_4326
    vertices = np.array(geom.exterior)
    # max x = east, max y = north
    east = max(vertices[:], key=lambda item:item[0])
    north = max(vertices[:], key=lambda item:item[1])
    west = min(vertices[:], key=lambda item:item[0])
    south = min(vertices[:], key=lambda item:item[1])
    east = east[0], east[1]
    north = north[0], north[1]
    west = west[0], west[1]
    south = south[0], south[1]
    df = df.append({'index': i, 'bid': row.bid, 'north': north, 'east': east, 'south': south, 'west': west}, ignore_index=True)

print("1")
# make new dataframes with only building id and point to concat into all sides df
gdf_north = df[['bid', 'north']].copy()
gdf_north['north'] = gdf_north['north'].apply(Point)
gdf_north['side'] = '360'
gdf_north = gpd.GeoDataFrame(gdf_north, geometry='north')
gdf_north = gdf_north.rename(columns={'north': 'point_3310'}).set_geometry('point_3310')
print("2")
gdf_south = df[['bid', 'south']].copy()
gdf_south['south'] = gdf_south['south'].apply(Point)
gdf_south['side'] = '180'
gdf_south = gpd.GeoDataFrame(gdf_south, geometry='south')
gdf_south = gdf_south.rename(columns={'south': 'point_3310'}).set_geometry('point_3310')
print("3")
gdf_east = df[['bid', 'east']].copy()
gdf_east['east'] = gdf_east['east'].apply(Point)
gdf_east['side'] = '90'
gdf_east = gpd.GeoDataFrame(gdf_east, geometry='east')
gdf_east = gdf_east.rename(columns={'east': 'point_3310'}).set_geometry('point_3310')
print("4")
gdf_west = df[['bid', 'west']].copy()
gdf_west['west'] = gdf_west['west'].apply(Point)
gdf_west['side'] = '270'
gdf_west = gpd.GeoDataFrame(gdf_west, geometry='west')
gdf_west = gdf_west.rename(columns={'west': 'point_3310'}).set_geometry('point_3310')

del df
print("5")
sides = pd.concat([gdf_north, gdf_east, gdf_south, gdf_west], ignore_index=True)
# # sides = sides.drop_duplicates(subset=['bid', 'side'], keep='first')
crs = {'init': 'epsg:3310'}
sides = gpd.GeoDataFrame(sides, crs=crs, geometry='point_3310')
sides['uid'] = sides['bid'].astype(int)+sides['side'].astype(int)/1000
sides.tail(1000)
t1 = time.time()
print('Total run time = ', t1 - t0)

sides['x'] = sides.geometry.x
sides['y'] = sides.geometry.y

# kdtree fcn for two lists: use nA for self, nB for all non-self building points
def ckd_nearest(gdA, gdB, bcol):
    nA= np.array(list(zip(gdA.geometry.x, gdA.geometry.y)) )
    nB = np.array(list(zip(gdB.geometry.x, gdA.geometry.y)) )
    btree = cKDTree(nB)
    dist, idx = btree.query(nA,k=1)
    df = pd.DataFrame.from_dict({'distance': dist.astype(int), 'bcol' : gdB.loc[idx, bcol].values })
    return df

# initialize output df
df = pd.DataFrame({'bid': [], 'uid': [], 'distance': [], 'nearest_idx': []})

# loop through kdtree for each point while excluding same building in search
all_sides = np.array(list(zip(sides.geometry.x, sides.geometry.y)))


for i, row in sides.iterrows():
#     t0 = time.time()
    queryable = pd.DataFrame({'x': [], 'y': [], 'uid': []})
    allsides = np.array(list(zip(sides.geometry.x, sides.geometry.y, sides.uid)))#, sides.side, #building[2] == row.bid, 

    for building in allsides:
        if (int(str(building[2]).split('.')[0]) != row.bid): 
            queryable = queryable.append({'x': building[0], 'y': building[1], 'uid': building[2]}, ignore_index=True)
    nA = np.array(pd.DataFrame({'x': [row.x], 'y': [row.y], 'uid': [row.uid]}))   
    nB = np.array(queryable)
    btree = cKDTree(nB)
    dist, idx = btree.query(nA, k=1)
    print(row.uid, dist, idx)
    row.distance = dist
    row.nearest_i = idx
    df = df.append({'bid': row.bid, 'uid': row.uid, 'distance': dist, 'nearest_idx': idx}, ignore_index=True)

#     print('Total run time = ', time.time() - t0)
    
# preview data
df.tail(500)

# save results
destination = '/home/kwaechte/data/comstock/'
name='bldg_offsets'
df.to_csv(str(destination)+str(name)+str('.csv'))
print("Data exported to ", destination, ".")
