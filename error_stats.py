# Copyright 2018 Samuel B. Powell, Washinton University in St. Louis
# 
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
from pylab import *
import scipy.stats as stats
from stats_utils import *
from plot_utils import *
import numpy.lib.recfunctions as rf

from kent_distribution import kent_mle
from mpl_toolkits.basemap import Basemap
from geographiclib.geodesic import Geodesic
WGS84 = Geodesic.WGS84

make_plots = False
make_maps = False
write_tables = False
#### Load Data ####

#we need to load:
# stats-nr-L1, stats-knn-2016-L1, stats-knn-2015-L1, stats-knn-2015-hi-L1

def load_csv(fname):
    #load csv into structured numpy array
    #assume first line is field names, autodetect dtypes for each column, remove whitespace
    x = genfromtxt(fname,delimiter=',',names=True,dtype=None,autostrip=True)
    #return as a record array to make access easier
    return x.view(np.recarray)

def write_csv(fname,rarry):
    delim = ', '
    newline = '\n'
    def _str(x):
        if isinstance(x,bytes): return str(x,'UTF-8')
        else: return str(x)
    with open(fname,'w') as f:
        #save recarray as csv file
        #header:
        f.write(delim.join(rarry.dtype.names))
        f.write(newline)
        #data:
        for row in rarry:
            f.write(delim.join(_str(c) for c in row))
            f.write(newline)


print('Loading data...',flush=True)
stats_naive = load_csv('../stats-nr-L1.csv')
stats_knn_16 = load_csv('../stats-knn-2016-L1.csv')
stats_knn_15 = load_csv('../stats-knn-2015-L1.csv')
stats_knn_hi_15 = load_csv('../stats-knn-2015-hi-L1.csv')

#the fields in each csv file are:
# name, site, time, lat, lon, depth, sun_head, sun_zen,
# sun_head_est, sun_zen_est, sun_pos_error, lat_est, lon_est, gps_error

#### Compute Errors ####
print('Computing errors...',flush=True)
#add heading and zenith error fields to each
stats_naive = rf.append_fields(stats_naive,['head_err','zen_err'],data=[[0.0]]*2,usemask=False,asrecarray=True)
stats_knn_16 = rf.append_fields(stats_knn_16,['head_err','zen_err'],data=[[0.0]]*2,usemask=False,asrecarray=True)
stats_knn_15 = rf.append_fields(stats_knn_15,['head_err','zen_err'],data=[[0.0]]*2,usemask=False,asrecarray=True)
stats_knn_hi_15 = rf.append_fields(stats_knn_hi_15,['head_err','zen_err'],data=[[0.0]]*2,usemask=False,asrecarray=True)

for s in (stats_naive, stats_knn_16, stats_knn_15, stats_knn_hi_15):
    #measured = truth + noise ==> noise = measured - truth
    s.head_err = angle_diff(s.sun_head_est, s.sun_head)
    s.zen_err = angle_diff(s.sun_zen_est, s.sun_zen)

#### Group Data ####
print('Grouping data...',flush=True)
stats_knn_hi_16 = stats_knn_16

#we want to compute for the naive and kNN GPS errors, heading errors, and zenith errors
#  RMS, paired T-stats, and p-values
#   grouped by site, by 2015 data, by 2016 data, and over all data
#   and with the same groupings for high-sun (sun > 40°)

#get naive stats for hi sun elevation (ie. low sun zenith)
stats_naive_hi = stats_naive[stats_naive.sun_zen < deg2rad(90-40)]

#get naive stats by sensor
stats_naive_15 = stats_naive[:len(stats_knn_15)]
stats_naive_16 = stats_naive[len(stats_knn_15):]

stats_naive_hi_15 = stats_naive_hi[:len(stats_knn_hi_15)]
stats_naive_hi_16 = stats_naive_hi[len(stats_knn_hi_15):]

#all knn data together:
stats_knn = rf.stack_arrays((stats_knn_15,stats_knn_16),autoconvert=True,asrecarray=True,usemask=False)
stats_knn_hi = rf.stack_arrays((stats_knn_hi_15,stats_knn_16),autoconvert=True,asrecarray=True,usemask=False)


#data from australia
sites_aus = [    
    b'Mermaid 15',
    b'Casuaria Dusk 15',
    b'No Name 15',
    b'Cobia Hole 15',
    b'Snake Pit 15',
    b'Casuaria Day 15',
    b'Mermaid (deep) 16',
    b'Mermaid 16',
    b'Fake Cobia Hole 16',
    b'Loomis 16',
]

aus_mask = sum([stats_knn['site'] == s for s in sites_aus],axis=0) > 0
aus_hi_mask = sum([stats_knn_hi['site'] == s for s in sites_aus],axis=0) > 0

stats_naive_aus = stats_naive[aus_mask]
stats_knn_aus = stats_knn[aus_mask]

stats_naive_hi_aus = stats_naive_hi[aus_hi_mask]
stats_knn_hi_aus = stats_knn_hi[aus_hi_mask]

#write csv files w/ errors
print('Writing stats files with errors...',flush=True)
write_csv('stats_naive.csv',stats_naive)
write_csv('stats_knn.csv',stats_knn)

#group by site, all sun elevations
sites = set(stats_knn.site)
stats_naive_by_site = {s:stats_naive[stats_naive.site == s] for s in sites}
stats_knn_by_site = {s:stats_knn[stats_knn.site == s] for s in sites}

#group by site, high sun elevation
sites_hi = set(stats_knn_hi.site)
stats_naive_hi_by_site = {s:stats_naive_hi[stats_naive_hi.site == s] for s in sites_hi}
stats_knn_hi_by_site = {s:stats_knn_hi[stats_knn_hi.site == s] for s in sites_hi}

#### Compute Statistics ####

#for kent distribution:
def gps2xyz(lat,lon):
    #convert latitude, longitude in degrees into xyz coordinates
    lat,lon = broadcast_arrays(deg2rad(lat),deg2rad(lon))
    xyz = empty(lat.shape+(3,))
    c = cos(lat)
    xyz[...,0] = c*cos(lon)
    xyz[...,1] = c*sin(lon)
    xyz[...,2] = sin(lat)
    return xyz

def xyz2gps(xyz,axis=-1):
    #convert x,y,z to latitude, longitude in degrees
    xyz = asarray(xyz)
    x,y,z = take(xyz,0,axis),take(xyz,1,axis),take(xyz,2,axis)
    lon = rad2deg(arctan2(y,x))
    r = hypot(x,y)
    lat = rad2deg(arctan2(z,r))
    return lat,lon

def gps_mean(lat,lon,axis=None):
    xyz = gps2xyz(lat,lon)
    return xyz2gps(mean(xyz.reshape((-1,3)),0))

def gps_dist(gps0,gps1):
    """gps_dist(gps0,gps1)
    Distance in m between gps coordinates
    gps0,gps1 : lat,lon pairs in degrees
    """
    return WGS84.Inverse(gps0[0],gps0[1],gps1[0],gps1[1])['s12']

def gps_centroid_dist(lats1, lons1, lats2, lons2):
    return gps_dist(gps_mean(lats1,lons2), gps_mean(lats2,lons2))

print('Computing stats...',flush=True)
#for all sun elevations
print('  for all sun elevations...',flush=True)
#compute RMS over various groupings:
gps_naive_rmse, head_naive_rmse, zen_naive_rmse = {}, {}, {}
gps_knn_rmse, head_knn_rmse, zen_knn_rmse = {}, {}, {}

gps_naive_me, head_naive_me = {}, {}
gps_knn_me, head_knn_me = {}, {}

for s in sites:
    sn = stats_naive_by_site[s]
    gps_naive_rmse[s] = rms(sn.gps_error)
    head_naive_rmse[s] = angle_rms(sn.head_err)
    zen_naive_rmse[s] = angle_rms(sn.zen_err)

    head_naive_me[s] = angle_mean(sn.head_err)
    gps_naive_me[s] = gps_centroid_dist(sn.lat, sn.lon, sn.lat_est, sn.lon_est)

    sk = stats_knn_by_site[s]
    gps_knn_rmse[s] = rms(sk.gps_error)
    head_knn_rmse[s] = angle_rms(sk.head_err)
    zen_knn_rmse[s] = angle_rms(sk.zen_err)

    head_knn_me[s] = angle_mean(sk.head_err)
    gps_knn_me[s] = gps_centroid_dist(sk.lat, sk.lon, sk.lat_est, sk.lon_est)

#by sensor
gps_naive_rmse[b'sensor 15'] = rms(stats_naive_15.gps_error)
head_naive_rmse[b'sensor 15'] = angle_rms(stats_naive_15.head_err)
zen_naive_rmse[b'sensor 15'] = angle_rms(stats_naive_15.zen_err)
head_naive_me[b'sensor 15'] = angle_mean(stats_naive_15.head_err)
gps_naive_me[b'sensor 15'] = nan #nan because this doesn't make sense

gps_knn_rmse[b'sensor 15'] = rms(stats_knn_15.gps_error)
head_knn_rmse[b'sensor 15'] = angle_rms(stats_knn_15.head_err)
zen_knn_rmse[b'sensor 15'] = angle_rms(stats_knn_15.zen_err)
head_knn_me[b'sensor 15'] = angle_mean(stats_knn_15.head_err)
gps_knn_me[b'sensor 15'] = nan

gps_naive_rmse[b'sensor 16'] = rms(stats_naive_16.gps_error)
head_naive_rmse[b'sensor 16'] = angle_rms(stats_naive_16.head_err)
zen_naive_rmse[b'sensor 16'] = angle_rms(stats_naive_16.zen_err)
head_naive_me[b'sensor 16'] = angle_mean(stats_naive_16.head_err)
gps_naive_me[b'sensor 16'] = nan

gps_knn_rmse[b'sensor 16'] = rms(stats_knn_16.gps_error)
head_knn_rmse[b'sensor 16'] = angle_rms(stats_knn_16.head_err)
zen_knn_rmse[b'sensor 16'] = angle_rms(stats_knn_16.zen_err)
head_knn_me[b'sensor 16'] = angle_mean(stats_knn_16.head_err)
gps_knn_me[b'sensor 16'] = nan

#australia, lizard
gps_naive_rmse[b'lizard island'] = rms(stats_naive_aus.gps_error)
head_naive_rmse[b'lizard island'] = angle_rms(stats_naive_aus.head_err)
zen_naive_rmse[b'lizard island'] = angle_rms(stats_naive_aus.zen_err)
head_naive_me[b'lizard island'] = angle_mean(stats_naive_aus.head_err)
gps_naive_me[b'lizard island'] = gps_centroid_dist(stats_naive_aus.lat, stats_naive_aus.lon, stats_naive_aus.lat_est, stats_naive_aus.lon_est)

gps_knn_rmse[b'lizard island'] = rms(stats_knn_aus.gps_error)
head_knn_rmse[b'lizard island'] = angle_rms(stats_knn_aus.head_err)
zen_knn_rmse[b'lizard island'] = angle_rms(stats_knn_aus.zen_err)
head_knn_me[b'lizard island'] = angle_mean(stats_knn_aus.head_err)
gps_knn_me[b'lizard island'] = gps_centroid_dist(stats_knn_aus.lat, stats_knn_aus.lon, stats_knn_aus.lat_est, stats_knn_aus.lon_est)


#all
gps_naive_rmse[b'all'] = rms(stats_naive.gps_error)
head_naive_rmse[b'all'] = angle_rms(stats_naive.head_err)
zen_naive_rmse[b'all'] = angle_rms(stats_naive.zen_err)
head_naive_me[b'all'] = angle_mean(stats_naive.head_err)
gps_naive_me[b'all'] = nan

gps_knn_rmse[b'all'] = rms(stats_knn.gps_error)
head_knn_rmse[b'all'] = angle_rms(stats_knn.head_err)
zen_knn_rmse[b'all'] = angle_rms(stats_knn.zen_err)
head_knn_me[b'all'] = angle_mean(stats_knn.head_err)
gps_knn_me[b'all'] = nan

#paired t-tests between naive & kNN
gps_tp, head_tp, zen_tp = {}, {}, {}

for s in sites:
    sn = stats_naive_by_site[s]
    sk = stats_knn_by_site[s]
    gps_tp[s] = stats.ttest_rel(sn.gps_error, sk.gps_error)
    head_tp[s] = stats.ttest_1samp(angle_diff(sn.head_err, sk.head_err),0)
    zen_tp[s] = stats.ttest_1samp(angle_diff(sn.zen_err, sk.zen_err),0)

#by sensor
gps_tp[b'sensor 15'] = stats.ttest_rel(stats_naive_15.gps_error, stats_knn_15.gps_error)
head_tp[b'sensor 15'] = stats.ttest_1samp(angle_diff(stats_naive_15.head_err, stats_knn_15.head_err),0)
zen_tp[b'sensor 15'] = stats.ttest_1samp(angle_diff(stats_naive_15.zen_err, stats_knn_15.zen_err),0)

gps_tp[b'sensor 16'] = stats.ttest_rel(stats_naive_16.gps_error, stats_knn_16.gps_error)
head_tp[b'sensor 16'] = stats.ttest_1samp(angle_diff(stats_naive_16.head_err, stats_knn_16.head_err),0)
zen_tp[b'sensor 16'] = stats.ttest_1samp(angle_diff(stats_naive_16.zen_err, stats_knn_16.zen_err),0)

#lizard
gps_tp[b'lizard island'] = stats.ttest_rel(stats_naive_aus.gps_error, stats_knn_aus.gps_error)
head_tp[b'lizard island'] = stats.ttest_1samp(angle_diff(stats_naive_aus.head_err, stats_knn_aus.head_err),0)
zen_tp[b'lizard island'] = stats.ttest_1samp(angle_diff(stats_naive_aus.zen_err, stats_knn_aus.zen_err),0)

#all
gps_tp[b'all'] = stats.ttest_rel(stats_naive.gps_error, stats_knn.gps_error)
head_tp[b'all'] = stats.ttest_1samp(angle_diff(stats_naive.head_err, stats_knn.head_err),0)
zen_tp[b'all'] = stats.ttest_1samp(angle_diff(stats_naive.zen_err, stats_knn.zen_err),0)

## for high sun elevations
print('  for high sun elevations...',flush=True)
#compute RMS over various groupings:
gps_naive_hi_rmse, head_naive_hi_rmse, zen_naive_hi_rmse = {}, {}, {}
gps_knn_hi_rmse, head_knn_hi_rmse, zen_knn_hi_rmse = {}, {}, {}

gps_naive_hi_me, head_naive_hi_me = {}, {}
gps_knn_hi_me, head_knn_hi_me = {}, {}

for s in sites_hi:
    sn = stats_naive_hi_by_site[s]
    gps_naive_hi_rmse[s] = rms(sn.gps_error)
    head_naive_hi_rmse[s] = angle_rms(sn.head_err)
    zen_naive_hi_rmse[s] = angle_rms(sn.zen_err)

    head_naive_hi_me[s] = angle_mean(sn.head_err)
    gps_naive_hi_me[s] = gps_centroid_dist(sn.lat, sn.lon, sn.lat_est, sn.lon_est)

    sk = stats_knn_hi_by_site[s]
    gps_knn_hi_rmse[s] = rms(sk.gps_error)
    head_knn_hi_rmse[s] = angle_rms(sk.head_err)
    zen_knn_hi_rmse[s] = angle_rms(sk.zen_err)

    head_knn_hi_me[s] = angle_mean(sk.head_err)
    gps_knn_hi_me[s] = gps_centroid_dist(sk.lat, sk.lon, sk.lat_est, sk.lon_est)

#by sensor
gps_naive_hi_rmse[b'sensor 15'] = rms(stats_naive_hi_15.gps_error)
head_naive_hi_rmse[b'sensor 15'] = angle_rms(stats_naive_hi_15.head_err)
zen_naive_hi_rmse[b'sensor 15'] = angle_rms(stats_naive_hi_15.zen_err)
head_naive_hi_me[b'sensor 15'] = angle_mean(stats_naive_hi_15.head_err)
gps_naive_hi_me[b'sensor 15'] = nan #nan because this doesn't make sense

gps_knn_hi_rmse[b'sensor 15'] = rms(stats_knn_hi_15.gps_error)
head_knn_hi_rmse[b'sensor 15'] = angle_rms(stats_knn_hi_15.head_err)
zen_knn_hi_rmse[b'sensor 15'] = angle_rms(stats_knn_hi_15.zen_err)
head_knn_hi_me[b'sensor 15'] = angle_mean(stats_knn_hi_15.head_err)
gps_knn_hi_me[b'sensor 15'] = nan

gps_naive_hi_rmse[b'sensor 16'] = rms(stats_naive_hi_16.gps_error)
head_naive_hi_rmse[b'sensor 16'] = angle_rms(stats_naive_hi_16.head_err)
zen_naive_hi_rmse[b'sensor 16'] = angle_rms(stats_naive_hi_16.zen_err)
head_naive_hi_me[b'sensor 16'] = angle_mean(stats_naive_hi_16.head_err)
gps_naive_hi_me[b'sensor 16'] = nan

gps_knn_hi_rmse[b'sensor 16'] = rms(stats_knn_hi_16.gps_error)
head_knn_hi_rmse[b'sensor 16'] = angle_rms(stats_knn_hi_16.head_err)
zen_knn_hi_rmse[b'sensor 16'] = angle_rms(stats_knn_hi_16.zen_err)
head_knn_hi_me[b'sensor 16'] = angle_mean(stats_knn_hi_16.head_err)
gps_knn_hi_me[b'sensor 16'] = nan

#lizard island
gps_naive_hi_rmse[b'lizard island'] = rms(stats_naive_hi_aus.gps_error)
head_naive_hi_rmse[b'lizard island'] = angle_rms(stats_naive_hi_aus.head_err)
zen_naive_hi_rmse[b'lizard island'] = angle_rms(stats_naive_hi_aus.zen_err)
head_naive_hi_me[b'lizard island'] = angle_mean(stats_naive_hi_aus.head_err)
gps_naive_hi_me[b'lizard island'] = gps_centroid_dist(stats_naive_hi_aus.lat, stats_naive_hi_aus.lon, stats_naive_hi_aus.lat_est, stats_naive_hi_aus.lon_est)

gps_knn_hi_rmse[b'lizard island'] = rms(stats_knn_hi_aus.gps_error)
head_knn_hi_rmse[b'lizard island'] = angle_rms(stats_knn_hi_aus.head_err)
zen_knn_hi_rmse[b'lizard island'] = angle_rms(stats_knn_hi_aus.zen_err)
head_knn_hi_me[b'lizard island'] = angle_mean(stats_knn_hi_aus.head_err)
gps_knn_hi_me[b'lizard island'] = gps_centroid_dist(stats_knn_hi_aus.lat, stats_knn_hi_aus.lon, stats_knn_hi_aus.lat_est, stats_knn_hi_aus.lon_est)

#all
gps_naive_hi_rmse[b'all'] = rms(stats_naive_hi.gps_error)
head_naive_hi_rmse[b'all'] = angle_rms(stats_naive_hi.head_err)
zen_naive_hi_rmse[b'all'] = angle_rms(stats_naive_hi.zen_err)
head_naive_hi_me[b'all'] = angle_mean(stats_naive_hi.head_err)
gps_naive_hi_me[b'all'] = nan

gps_knn_hi_rmse[b'all'] = rms(stats_knn_hi.gps_error)
head_knn_hi_rmse[b'all'] = angle_rms(stats_knn_hi.head_err)
zen_knn_hi_rmse[b'all'] = angle_rms(stats_knn_hi.zen_err)
head_knn_hi_me[b'all'] = angle_mean(stats_knn_hi.head_err)
gps_knn_hi_me[b'all'] = nan

#paired t-tests between naive_hi & knn_hi
gps_hi_tp, head_hi_tp, zen_hi_tp = {}, {}, {}

for s in sites_hi:
    sn = stats_naive_hi_by_site[s]
    sk = stats_knn_hi_by_site[s]
    gps_hi_tp[s] = stats.ttest_rel(sn.gps_error, sk.gps_error)
    head_hi_tp[s] = stats.ttest_1samp(angle_diff(sn.head_err, sk.head_err),0)
    zen_hi_tp[s] = stats.ttest_1samp(angle_diff(sn.zen_err, sk.zen_err),0)

#by sensor
gps_hi_tp[b'sensor 15'] = stats.ttest_rel(stats_naive_hi_15.gps_error, stats_knn_hi_15.gps_error)
head_hi_tp[b'sensor 15'] = stats.ttest_1samp(angle_diff(stats_naive_hi_15.head_err, stats_knn_hi_15.head_err),0)
zen_hi_tp[b'sensor 15'] = stats.ttest_1samp(angle_diff(stats_naive_hi_15.zen_err, stats_knn_hi_15.zen_err),0)

gps_hi_tp[b'sensor 16'] = stats.ttest_rel(stats_naive_hi_16.gps_error, stats_knn_hi_16.gps_error)
head_hi_tp[b'sensor 16'] = stats.ttest_1samp(angle_diff(stats_naive_hi_16.head_err, stats_knn_hi_16.head_err),0)
zen_hi_tp[b'sensor 16'] = stats.ttest_1samp(angle_diff(stats_naive_hi_16.zen_err, stats_knn_hi_16.zen_err),0)

#lizard
gps_hi_tp[b'lizard island'] = stats.ttest_rel(stats_naive_hi_aus.gps_error, stats_knn_hi_aus.gps_error)
head_hi_tp[b'lizard island'] = stats.ttest_1samp(angle_diff(stats_naive_hi_aus.head_err, stats_knn_hi_aus.head_err),0)
zen_hi_tp[b'lizard island'] = stats.ttest_1samp(angle_diff(stats_naive_hi_aus.zen_err, stats_knn_hi_aus.zen_err),0)

#all
gps_hi_tp[b'all'] = stats.ttest_rel(stats_naive_hi.gps_error, stats_knn_hi.gps_error)
head_hi_tp[b'all'] = stats.ttest_1samp(angle_diff(stats_naive_hi.head_err, stats_knn_hi.head_err),0)
zen_hi_tp[b'all'] = stats.ttest_1samp(angle_diff(stats_naive_hi.zen_err, stats_knn_hi.zen_err),0)

#### Generate Tables and Charts ####
print('Making tables and figures...',flush=True)
#better site names:
display_sites = {
    b'Mermaid 15': 'Mermaid 1',
    b'Casuaria Dusk 15': 'Casuaria 1',
    b'No Name 15': 'No Name',
    b'Cobia Hole 15': 'Cobia Hole 1',
    b'Snake Pit 15': 'Snake Pit',
    b'Casuaria Day 15': 'Casuaria 2',
    b'Finland 15': 'Tvärminne',
    b'Electric Beach 15': 'Electric Beach',
    b'Miami 16': 'Miami',
    b'Mermaid (deep) 16': 'Mermaid 2',
    b'Mermaid 16': 'Mermaid 3',
    b'Fake Cobia Hole 16': 'Cobia Hole 2',
    b'Loomis 16': 'Loomis',
    b'lizard island': 'Lizard Island',
    b'sensor 15': 'All Sensor 1',
    b'sensor 16': 'All Sensor 2',
    b'all': 'All',
}

#chronological order:
site_order = [    
    b'Mermaid 15',
    b'Casuaria Dusk 15',
    b'No Name 15',
    b'Cobia Hole 15',
    b'Snake Pit 15',
    b'Casuaria Day 15',
    b'Finland 15',
    b'Electric Beach 15',
    b'Miami 16',
    b'Mermaid (deep) 16',
    b'Mermaid 16',
    b'Fake Cobia Hole 16',
    b'Loomis 16',
    b'lizard island',
    b'sensor 15',
    b'sensor 16',
    b'all',
]

site_order_limited = [
    b'Finland 15',
    b'Electric Beach 15',
    b'Miami 16',
    b'lizard island',
    #b'sensor 15',
    #b'sensor 16',
    b'all',
]
sol_idx = [site_order.index(s) for s in site_order_limited]

site_order_hi = [s for s in site_order if s in sites_hi]
site_order_hi += [b'lizard island',b'sensor 15', b'sensor 16', b'all']

site_order_limited_hi = [s for s in site_order_limited if s in sites_hi]
site_order_limited_hi += [b'lizard island',b'all']
sol_idx_hi = [site_order_hi.index(s) for s in site_order_limited_hi]

print('  all data...',flush=True)
print('    rms error table...',flush=True)
headers = [
    'Site',
    'Naive GPS RMSE (km)',
    'Naive Heading RMSE (deg)',
    'Naive Elevation RMSE (deg)',
    'Naive GPS ME (km)',
    'Naive Heading ME (deg)',
    'kNN GPS RMSE (km)',
    'kNN Heading Error (deg)',
    'kNN Elevation Error (deg)',
    'kNN GPS ME (km)',
    'kNN Heading ME (deg)',
    'GPS p-value',
    'Heading p-value',
    'Elevation p-value'
]
rms_table = [headers]
for site in site_order:
    row = [display_sites[site]]
    row += [gps_naive_rmse[site]/1000, rad2deg(head_naive_rmse[site])]
    row += [rad2deg(zen_naive_rmse[site]), gps_naive_me[site]/1000]
    row += [rad2deg(head_naive_me[site])]
    row += [gps_knn_rmse[site]/1000, rad2deg(head_knn_rmse[site])]
    row += [rad2deg(zen_knn_rmse[site]),gps_knn_me[site]/1000]
    row += [rad2deg(head_knn_me[site])]
    row += [gps_tp[site].pvalue,head_tp[site].pvalue,zen_tp[site].pvalue]
    rms_table.append(row)

if write_tables:
    with open('summary_errors_by_site.csv','w') as f:
        for row in rms_table:
            f.write(', '.join(str(c) for c in row))
            f.write('\n')

def barchart(ax, sites, nv, knn, pvals):
    n = len(sites)
    bar_width = 0.4
    bar_margin = 0
    color_a = 'r'
    color_b = 'b'
    xx = arange(n)+1
    xa = xx - bar_width - bar_margin/2
    xb = xa + bar_width + bar_margin

    bars_a = ax.bar(xa,nv,bar_width,0,color=color_a)
    bars_b = ax.bar(xb,knn,bar_width,0,color=color_b)

    yrange = ax.get_ylim()[1]
    texts = []
    #place a '*' above sites w/ p < 0.05
    for x,a,b,p in zip(xx,nv,knn,pvals):
        if p < 0.05:
            h = max(0, a, b)+yrange*0.01
            texts.append(ax.text(x,h,'*',ha='center',va='bottom'))
    
    #what is the maximum y of the text bounding boxes?
    if texts:
        ymax = max(text_bbox(t).ymax for t in texts)
        text_height = text_bbox(texts[0]).height
        if ymax > yrange:
            ax.set_ylim(0,ymax + text_height)

    ax.set_xticks(arange(n)+1)
    ax.set_xlim(0,n+1)
    ax.set_xticklabels(sites,ha='right',rotation=25)
    ax.legend((bars_a[0],bars_b[0]),('Naive','kNN'))

    return bars_a, bars_b

labels = [display_sites[s] for s in site_order]
labels_limited = [display_sites[s] for s in site_order_limited]

gps_n_rmse = [gps_naive_rmse[s]/1000 for s in site_order]
gps_k_rmse = [gps_knn_rmse[s]/1000 for s in site_order]
gps_n_me = [gps_naive_me[s]/1000 for s in site_order]
gps_k_me = [gps_knn_me[s]/1000 for s in site_order]
gps_p = [gps_tp[s].pvalue for s in site_order]

head_n_rmse = [rad2deg(head_naive_rmse[s]) for s in site_order]
head_k_rmse = [rad2deg(head_knn_rmse[s]) for s in site_order]
head_n_me = [rad2deg(head_naive_me[s]) for s in site_order]
head_k_me = [rad2deg(head_knn_me[s]) for s in site_order]
head_p = [head_tp[s].pvalue for s in site_order]

zen_n_rmse = [rad2deg(zen_naive_rmse[s]) for s in site_order]
zen_k_rmse = [rad2deg(zen_knn_rmse[s]) for s in site_order]
zen_p = [zen_tp[s].pvalue for s in site_order]

if make_plots:
    print('    bar charts...',flush=True)
    #fig,ax = subplots(num='GPS RMSE all')
    #barchart(ax,labels,gps_n_rmse,gps_k_rmse,gps_p)
    #ax.set_ylabel('RMS GPS Error (km)')
    #fig.subplots_adjust(0.1,0.13,0.9,0.98)

    fig,ax = subplots(num='GPS RMSE')
    barchart(ax,labels_limited,take(gps_n_rmse,sol_idx),take(gps_k_rmse,sol_idx),take(gps_p,sol_idx))
    ax.set_ylabel('RMS GPS Error (km)')
    fig.subplots_adjust(0.1,0.13,0.9,0.98)

    #fig,ax = subplots(num='GPS ME all')
    #barchart(ax,labels[:-3],gps_n_me[:-3],gps_k_me[:-3],gps_p[:-3])
    #ax.set_ylabel('GPS Error (km)')
    #fig.subplots_adjust(0.1,0.13,0.9,0.98)

    fig,ax = subplots(num='GPS ME')
    barchart(ax,labels_limited[:-1],take(gps_n_me,sol_idx[:-1]),take(gps_k_me,sol_idx[:-1]),take(gps_p,sol_idx[:-1]))
    ax.set_ylabel('GPS Error (km)')
    fig.subplots_adjust(0.1,0.13,0.9,0.98)

    #fig,ax = subplots(num='Heading RMSE all')
    #barchart(ax,labels,head_n_rmse,head_k_rmse,head_p)
    #myticks(labels='{:.0f}°',axes=ax)
    #ax.set_ylabel('RMS Heading Error')
    #fig.subplots_adjust(0.1,0.13,0.9,0.98)

    fig,ax = subplots(num='Heading RMSE')
    barchart(ax,labels_limited,take(head_n_rmse,sol_idx),take(head_k_rmse,sol_idx),take(head_p,sol_idx))
    myticks(labels='{:.0f}°',axes=ax)
    ax.set_ylabel('RMS Heading Error')
    fig.subplots_adjust(0.1,0.13,0.9,0.98)

    #fig,ax = subplots(num='Heading ME all')
    #barchart(ax,labels,head_n_me,head_k_me,head_p)
    #myticks(labels='{:.0f}°',axes=ax)
    #ax.set_ylabel('Heading Error')
    #fig.subplots_adjust(0.1,0.13,0.9,0.98)

    fig,ax = subplots(num='Heading ME')
    barchart(ax,labels_limited,take(head_n_me,sol_idx),take(head_k_me,sol_idx),take(head_p,sol_idx))
    myticks(labels='{:.0f}°',axes=ax)
    ax.set_ylabel('Heading Error')
    fig.subplots_adjust(0.1,0.13,0.9,0.98)

    #fig,ax = subplots(num='Elevation RMSE all')
    #barchart(ax,labels,zen_n_rmse,zen_k_rmse,zen_p)
    #myticks(labels='{:.0f}°',axes=ax)
    #ax.set_ylabel('RMS Elevation Error')
    #fig.subplots_adjust(0.1,0.13,0.9,0.98)

    fig,ax = subplots(num='Elevation RMSE')
    barchart(ax,labels_limited,take(zen_n_rmse,sol_idx),take(zen_k_rmse,sol_idx),take(zen_p,sol_idx))
    myticks(labels='{:.0f}°',axes=ax)
    ax.set_ylabel('RMS Elevation Error')
    fig.subplots_adjust(0.1,0.13,0.9,0.98)

#high sun
print('  high sun...',flush=True)
print('    rms error table...',flush=True)
rms_table_hi = [headers]
for site in site_order_hi:
    row  = [display_sites[site]]
    row += [gps_naive_hi_rmse[site]/1000, rad2deg(head_naive_hi_rmse[site])]
    row += [rad2deg(zen_naive_hi_rmse[site]), gps_naive_hi_me[site]/1000]
    row += [rad2deg(head_naive_hi_me[site])]
    row += [gps_knn_hi_rmse[site]/1000, rad2deg(head_knn_hi_rmse[site])]
    row += [rad2deg(zen_knn_hi_rmse[site]),gps_knn_hi_me[site]/1000]
    row += [rad2deg(head_knn_hi_me[site])]
    row += [gps_hi_tp[site].pvalue,head_hi_tp[site].pvalue,zen_hi_tp[site].pvalue]
    rms_table_hi.append(row)

if write_tables:
    with open('summary_errors_by_site_high_sun.csv','w') as f:
        for row in rms_table_hi:
            f.write(', '.join(str(c) for c in row))
            f.write('\n')

labels_hi = [display_sites[s] for s in site_order_hi]
labels_limited_hi = [display_sites[s] for s in site_order_limited_hi]

gps_hi_n_rmse = [gps_naive_hi_rmse[s]/1000 for s in site_order_hi]
gps_hi_k_rmse = [gps_knn_hi_rmse[s]/1000 for s in site_order_hi]
gps_hi_n_me = [gps_naive_hi_me[s]/1000 for s in site_order_hi]
gps_hi_k_me = [gps_knn_hi_me[s]/1000 for s in site_order_hi]
gps_hi_p = [gps_hi_tp[s].pvalue for s in site_order_hi]

head_hi_n_rmse = [rad2deg(head_naive_hi_rmse[s]) for s in site_order_hi]
head_hi_k_rmse = [rad2deg(head_knn_hi_rmse[s]) for s in site_order_hi]
head_hi_n_me = [rad2deg(head_naive_hi_me[s]) for s in site_order_hi]
head_hi_k_me = [rad2deg(head_knn_hi_me[s]) for s in site_order_hi]
head_hi_p = [head_hi_tp[s].pvalue for s in site_order_hi]

zen_hi_n_rmse = [rad2deg(zen_naive_hi_rmse[s]) for s in site_order_hi]
zen_hi_k_rmse = [rad2deg(zen_knn_hi_rmse[s]) for s in site_order_hi]
zen_hi_p = [zen_hi_tp[s].pvalue for s in site_order_hi]

if make_plots:
    print('    bar charts...',flush=True)
    #fig,ax = subplots(num='GPS RMSE all (High Sun)')
    #barchart(ax,labels_hi,gps_hi_n_rmse,gps_hi_k_rmse,gps_hi_p)
    #ax.set_ylabel('RMS GPS Error (km)')
    #fig.subplots_adjust(0.1,0.13,0.9,0.98)

    #fig,ax = subplots(num='GPS ME all (High Sun)')
    #barchart(ax,labels_hi[:-3],gps_hi_n_me[:-3],gps_hi_k_me[:-3],gps_hi_p[:-3])
    #ax.set_ylabel('GPS Error (km)')
    #fig.subplots_adjust(0.1,0.13,0.9,0.98)

    #fig,ax = subplots(num='Heading RMSE all (High Sun)')
    #barchart(ax,labels_hi,head_hi_n_rmse,head_hi_k_rmse,head_hi_p)
    #ax.set_ylabel('RMS Heading Error')
    #myticks(labels='{:.0f}°',axes=ax)
    #fig.subplots_adjust(0.1,0.13,0.9,0.98)

    #fig,ax = subplots(num='Heading ME all (High Sun)')
    #barchart(ax,labels_hi,head_hi_n_me,head_hi_k_me,head_hi_p)
    #ax.set_ylabel('Heading Error')
    #myticks(labels='{:.0f}°',axes=ax)
    #fig.subplots_adjust(0.1,0.13,0.9,0.98)

    #fig,ax = subplots(num='Elevation Error all (High Sun)')
    #barchart(ax,labels_hi,zen_hi_n_rmse,zen_hi_k_rmse,zen_hi_p)
    #ax.set_ylabel('RMS Elevation Error')
    #myticks(labels='{:.0f}°',axes=ax)
    #fig.subplots_adjust(0.1,0.13,0.9,0.98)

    fig,ax = subplots(num='GPS RMSE (High Sun)')
    barchart(ax,labels_limited_hi,take(gps_hi_n_rmse,sol_idx_hi),take(gps_hi_k_rmse,sol_idx_hi),take(gps_hi_p,sol_idx_hi))
    ax.set_ylabel('RMS GPS Error (km)')
    fig.subplots_adjust(0.1,0.13,0.9,0.98)

    fig,ax = subplots(num='GPS ME (High Sun)')
    barchart(ax,labels_limited_hi[:-1],take(gps_hi_n_me,sol_idx_hi[:-1]),take(gps_hi_k_me,sol_idx_hi[:-1]),take(gps_hi_p,sol_idx_hi[:-1]))
    ax.set_ylabel('GPS Error (km)')
    fig.subplots_adjust(0.1,0.13,0.9,0.98)

    fig,ax = subplots(num='Heading RMSE (High Sun)')
    barchart(ax,labels_limited_hi,take(head_hi_n_rmse,sol_idx_hi),take(head_hi_k_rmse,sol_idx_hi),take(head_hi_p,sol_idx_hi))
    ax.set_ylabel('RMS Heading Error')
    myticks(labels='{:.0f}°',axes=ax)
    fig.subplots_adjust(0.1,0.13,0.9,0.98)

    fig,ax = subplots(num='Heading ME (High Sun)')
    barchart(ax,labels_limited_hi,take(head_hi_n_me,sol_idx_hi),take(head_hi_k_me,sol_idx_hi),take(head_hi_p,sol_idx_hi))
    ax.set_ylabel('Heading Error')
    myticks(labels='{:.0f}°',axes=ax)
    fig.subplots_adjust(0.1,0.13,0.9,0.98)

    fig,ax = subplots(num='Elevation Error (High Sun)')
    barchart(ax,labels_limited_hi,take(zen_hi_n_rmse,sol_idx_hi),take(zen_hi_k_rmse,sol_idx_hi),take(zen_hi_p,sol_idx_hi))
    ax.set_ylabel('RMS Elevation Error')
    myticks(labels='{:.0f}°',axes=ax)
    fig.subplots_adjust(0.1,0.13,0.9,0.98)


## Maps ##
if make_maps:
    print('Generating Maps...',flush=True)
#all australia, and per site maps:
#show true location, naive estimate, knn estimate



def map_line(lat1,lon1,lat2,lon2,n=10):
    #solve geodesic between positions
    gd = WGS84.Inverse(lat1,lon1,lat2,lon2)
    #create line from pos 1 to 2
    line = WGS84.Line(lat1,lon1,gd['azi1'])
    #evaluate line at n points
    arcdists = gd['a12']*linspace(0,1,n)[1:-1]
    lats,lons = [lat1],[lon1]
    for ad in arcdists:
        p = line.ArcPosition(ad)
        lats.append(p['lat2'])
        lons.append(p['lon2'])
    lats.append(lat2)
    lons.append(lon2)
    return array(lats), array(lons)

def plot_map_vs(ax, size, stats_naive, stats_knn, res='c'):
    
    #true lat, lon
    lat = angle_mean(stats_knn['lat'],period=360)
    lon = angle_mean(stats_knn['lon'],period=360)
    
    w, h = size
    lat_min, lat_max = max(lat-h/2,-85), min(lat+h/2,85)
    lon_min, lon_max = lon-w/2, lon+w/2
    mp = Basemap(lon_min,lat_min,lon_max,lat_max,resolution=res,ax=ax)
    #gray continents, blue water & lakes
    mp.drawmapboundary(fill_color='#bbf2f2')
    mp.fillcontinents(color='#cccccc',lake_color='#bbf2f2')

    #draw darker blue grid under the continent
    mp.drawparallels(arange(lat_min,lat_max+5,10).round(-1),dashes=[10,10],labels=[1,1,1,1],zorder=0,color='#99cccc')
    mp.drawmeridians(angle_wrap(arange(lon_min,lon_max+5,10).round(-1),360),dashes=[10,10],labels=[1,1,1,1],zorder=0,color='#99cccc')

    #draw lines from knn to naive
    for sn, sk in zip(stats_naive,stats_knn):
        latn, lonn = sn['lat_est'], sn['lon_est']
        latk, lonk = sk['lat_est'], sk['lon_est']
        lt,ln = map_line(latn,lonn,latk,lonk)
        mp.plot(ln,lt,'k:',lw=0.5,latlon=True)
        mp.plot(lonn,latn,'rs',latlon=True)
        mp.plot(lonk,latk,'b.',latlon=True)
    #draw true location last
    mp.plot(lon,lat,'gD',latlon=True)

def plot_on_map(mp, stats_map, kent=True, plot_points=True,wrap_lon=False,scale=1,truth_style=None,mean_style=None,est_style=None,pdf_style=None):

    if truth_style is None:
        truth_style = {'marker':'+','mfc':'None','mec':'b','mew':2*scale,'ms':7*scale}
    if mean_style is None:
        mean_style = {'marker':'x','mfc':'None','mec':'r','mew':2*scale,'ms':7*scale}
    if est_style is None:
        est_style = {'ls':'None','marker':'.','mfc':'k','mec':'None','ms':3*scale}
    if pdf_style is None:
        pdf_style = {'colors':'k','linewidths':0.25*scale}

    lat = angle_mean(stats_map['lat'],period=360)
    lon = angle_mean(stats_map['lon'],period=360)

    xyz = gps2xyz(stats_map['lat_est'],stats_map['lon_est'])
    if(len(stats_map) < 6): kent=False
    if kent:
        #use kent distribution to get mean and spread
        kd = kent_mle(xyz)
        #mean point:
        lat_m, lon_m = xyz2gps(kd.gamma1)
        #evaluate pdf over map so we can get contours
        w = abs(mp.lonmax-mp.lonmin)
        h = abs(mp.latmax-mp.latmin)
        lons,lats = mp.makegrid(round(w)+1,round(h)+1)
        pdf_g = kd.pdf(gps2xyz(lats,lons))
        #pdf levels corresponding to 3,2,1 stds
        #note that this is an approximation based on the Kent distribution being similar to the Normal distribution
        #if we didn't want to use this, we could numerically integrate the grid of pdf values to find the level sets
        # e.g: the solid angle of each grid cell is sr_g = cos(deg2rad(lats))*deg2rad(lats[0,0]-lats[1,0])*deg2rad(lons[0,0]-lons[0,1])
        #       sum all of the cells over a level: cmf(p) = sum((sr_g*pdf_g)[pdf_g>=p])
        # these seem to be within 1% or so of the expected integrated probability mass:
        pdf_levels = kd.pdf_max()*exp(-0.5*array([3,2,1])**2)
    else:
        #take mean, but don't try to estimate the distribution
        lat_m, lon_m = xyz2gps(mean(xyz,axis=0))
        pdf_g = None

    #draw estimate points
    if plot_points:
        mp.plot(stats_map['lon_est'], stats_map['lat_est'],latlon=True,**est_style)

    if wrap_lon:
        #draw true location
        mp.plot(wrap(lon,mp.lonmin,mp.lonmax),lat,latlon=True,**truth_style)
        #draw mean location
        mp.plot(wrap(lon_m,mp.lonmin,mp.lonmax),lat_m,latlon=True,**mean_style)
    else:
        #draw true location
        mp.plot(lon,lat,latlon=True,**truth_style)
        #draw mean location
        mp.plot(lon_m,lat_m,latlon=True,**mean_style)
    #draw contours for pdf
    if pdf_g is not None:
        mp.contour(lons,lats,pdf_g,pdf_levels,latlon=True,**pdf_style)


def plot_map(ax, size, stats_map, res='c', kent=True,**kw):
    #true lat, lon
    lat = angle_mean(stats_map['lat'],period=360)
    lon = angle_mean(stats_map['lon'],period=360)
    
    w, h = size
    lat_min, lat_max = max(lat-h/2,-85), min(lat+h/2,85)
    lon_min, lon_max = lon-w/2, lon+w/2
    mp = Basemap(lon_min,lat_min,lon_max,lat_max,resolution=res,ax=ax)
    #gray continents, blue water & lakes
    mp.drawmapboundary(fill_color='#bbf2f2')
    mp.fillcontinents(color='#cccccc',lake_color='#bbf2f2')

    #draw darker blue grid under the continent
    lineprops = {
        'labels':[1,1,1,1],
        'labelstyle':'+/-',
        'dashes':(None,None),
        'zorder':0,
        'color':'#99cccc',
        'linewidth':0.2,
    }
    mp.drawparallels(arange(lat_min,lat_max+5,10).round(-1),**lineprops)
    mp.drawmeridians(angle_wrap(arange(lon_min,lon_max+5,10).round(-1),360),**lineprops)

    plot_on_map(mp, stats_map, kent,**kw)
    return mp
 

if make_maps:

    print(' Australia, naive...',flush=True)
    fig,ax = subplots(num='map_naive_aus')
    plot_map(ax,(40,40),stats_naive_aus,'i') #i for final figure

    print(' Australia, kNN...',flush=True)
    fig,ax = subplots(num='map_knn_aus')
    plot_map(ax,(40,40),stats_knn_aus,'i') #i for final figure

    print(' Australia, high sun, naive...',flush=True)
    fig,ax = subplots(num='map_naive_hi_aus')
    plot_map(ax,(40,40),stats_naive_hi_aus,'i') #i for final figure

    print(' Australia, high sun, kNN...',flush=True)
    fig,ax = subplots(num='map_knn_hi_aus')
    plot_map(ax,(40,40),stats_knn_hi_aus,'i') #i for final figure

    for s in sites:
        ds = display_sites[s]
        print(' {}, naive...'.format(ds),flush=True)
        fig,ax = subplots(num='map_naive_{}'.format(display_sites[s]))
        plot_map(ax,(40,40),stats_naive_by_site[s],'i',False)

        print(' {}, kNN...'.format(ds),flush=True)
        fig,ax = subplots(num='map_knn_{}'.format(display_sites[s]))
        plot_map(ax,(40,40),stats_knn_by_site[s],'i',False)

    for s in sites_hi:
        print(' {}, high sun, naive...'.format(ds),flush=True)
        fig,ax = subplots(num='map_naive_hi_{}'.format(display_sites[s]))
        plot_map(ax,(40,40),stats_naive_hi_by_site[s],'i',False)

        print(' {}, high sun, kNN...'.format(ds),flush=True)
        fig,ax = subplots(num='map_knn_hi_{}'.format(display_sites[s]))
        plot_map(ax,(40,40),stats_knn_hi_by_site[s],'i',False)

    #global map
    fig,ax = subplots(num='map_naive_global',figsize=(13,6))
    mp = Basemap(-270,-85,90,85,resolution='c',ax=ax)
    mp.drawmapboundary(fill_color='#bbf2f2')
    mp.fillcontinents(color='#cccccc',lake_color='#bbf2f2')
    mp.drawparallels(arange(-180,181,15),dashes=(None,None),labels=[1,1,1,1],zorder=0,color='#99cccc',labelstyle='+/-')
    mp.drawmeridians(angle_wrap(arange(-180,181,30),360),dashes=(None,None),labels=[1,1,1,1],zorder=0,color='#99cccc',labelstyle='+/-');
    plot_on_map(mp,stats_naive_by_site[b'Finland 15'],kent=False,plot_points=True,wrap_lon=True)
    plot_on_map(mp,stats_naive_by_site[b'Miami 16'],kent=True,plot_points=True,wrap_lon=True)
    plot_on_map(mp,stats_naive_by_site[b'Electric Beach 15'],kent=True,plot_points=True,wrap_lon=True)
    plot_on_map(mp,stats_naive_aus,kent=True,plot_points=True,wrap_lon=True)
    tight_layout(rect=(0,0.05,1,0.95))

    fig,ax = subplots(num='map_naive_global',figsize=(13,6))
    mp = Basemap(-270,-85,90,85,resolution='c',ax=ax)
    mp.drawmapboundary(fill_color='#bbf2f2')
    mp.fillcontinents(color='#cccccc',lake_color='#bbf2f2')
    mp.drawparallels(arange(-180,181,15),dashes=(None,None),labels=[1,1,1,1],zorder=0,color='#99cccc',labelstyle='+/-')
    mp.drawmeridians(angle_wrap(arange(-180,181,30),360),dashes=(None,None),labels=[1,1,1,1],zorder=0,color='#99cccc',labelstyle='+/-');
    plot_on_map(mp,stats_naive_by_site[b'Finland 15'],kent=False,plot_points=True,wrap_lon=True)
    plot_on_map(mp,stats_naive_by_site[b'Miami 16'],kent=True,plot_points=True,wrap_lon=True)
    plot_on_map(mp,stats_naive_by_site[b'Electric Beach 15'],kent=True,plot_points=True,wrap_lon=True)
    plot_on_map(mp,stats_naive_aus,kent=True,plot_points=True,wrap_lon=True)
    tight_layout(rect=(0,0.05,1,0.95))





