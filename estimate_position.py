# Copyright 2018 Samuel B. Powell, Washington University in St. Louis
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
import h5py, os, sys, traceback
from glob import glob
from IPython.core.debugger import Tracer; set_trace = Tracer()

#my code:

#Bug
from polprocess import file_newer, load_cal, process

#TODO: sunpos_az empty? deleted arcdist
#from sunpos.csunpos import sunpos_az as sunpos

##TODO: rewrite import
#sunpos_az: Compute the azimuth and zenith angles of the sun as viewed at the given time and location.
#sunpos: Compute the observed and topocentric coordinates of the sun as viewed at the given time and location.
#sunpos_az is lieterally the cut version of sunpos with exactly the same input and output[...,:2]
from sunpos.csunpos import sunpos_az, sunpos
from simulate import oceanstokes, oceanaop
from polarization import *
from plot_utils import *
from stats_utils import *
from geomag.emm import EMMMesh

from datetime import datetime, timedelta, timezone
import scipy.stats as stats
from scipy.optimize import minimize
from scipy.spatial import cKDTree as KDTree
from geographiclib.geodesic import Geodesic
os.environ['PROJ_LIB'] = 'D:/NoSpace/Anaconda3/share/proj'
from mpl_toolkits.basemap import Basemap
from pint import UnitRegistry
import yaml
units = UnitRegistry()
WGS84 = Geodesic.WGS84

utc = timezone.utc
utcp10 = timezone(timedelta(hours=10))

#
def _p(v,l,*args,**kw):
    """print if v, prefix l centered dots"""
    if v:
        if l: print('\xB7'*l,end='') #B7 is centered dot
        print(*args,**kw)

class S2KDTree(KDTree):
    """Spatial data structure in S2 space, parameterized by azimuth, zenith (polar) angles"""
    def __init__(self,data_S2,leafsize=10):
        self._data_S2 = data_S2
        self._data_R3 = S2KDTree.S2toR3(data_S2)
        super().__init__(self._data_R3,leafsize)

    #TODO: added weights and cumulative
    def count_neighbors(self,other,r,p=2.0, weights=None, cumulative=None):
        """As KDTree.count_neighbors, except:
            other is an S2KDTree
            r is an arcdistance in radians
        """
        return super().count_neighbors(other,S2KDTree.AtoC(r),p)

    #TODO: added n_jobs
    def query(self,x,k=1,eps=0,p=2,distance_upper_bound=inf, n_jobs=1):
        """As KDTree.query, but all distances are arcdistance in radians"""
        if distance_upper_bound != inf:
            distance_upper_bound = S2KDTree.AtoC(distance_upper_bound)
        x = S2KDTree.S2toR3(x)
        d,i = super().query(x,k,eps,p,distance_upper_bound)
        return S2KDTree.CtoA(d),i
    def query_ball_point(self,x,r,p=2.0,eps=0):
        """As KDTree.query_ball_point, but r is an arcdistance in radians"""
        r = S2KDTree.AtoC(r)
        x = S2KDTree.S2toR3(x)
        return super().query_ball_point(x,r,p,eps)
    def query_ball_tree(self,other,r,p=2.0,eps=0):
        """As KDTree.query_ball_tree, except:
            other is an S2KDTree
            r is an arcdistance in radians
        """
        r = S2KDTree.AtoC(r)
        return super().query_ball_tree(other,r,p,eps)
    def query_pairs(self,r,p=2.0,eps=0):
        """As KDTree.query_pairs, but r is an arcdistance in radians"""
        r = S2KDTree.AtoC(r)
        return super().query_pairs(r,p,eps)
    def sparse_distance_matrix(self,other,max_distance,p=2.0):
        """As KDTree.sparse_distance_matrix, except:
            other is an S2KDTree
            max_distance is an arcdistance in radians
            the returned matrix is of arcdistances in radians
            """
        max_distance = S2KDTree.AtoC(max_distance)
        d = super().sparse_distance_matrix(other,max_distance,p)
        return S2KDTree.CtoA(d)
    @staticmethod
    def S2toR3(data_S2):
        """Transform S2 azimuth,zenith angle pairs to R3 points on the unit sphere
            x, y, z = cos(azimuth)sin(zenith), sin(azimuth)sin(zenith), cos(zenith)
        """
        data_S2 = atleast_2d(data_S2)
        if data_S2.shape[-1]!= 2:
            raise ValueError("data_S2 must have shape (...,2)")
        a,z = data_S2[...,0],data_S2[...,1]
        data_R3 = empty(data_S2.shape[:-1]+(3,))
        sz = sin(z)
        data_R3[...,0] = cos(a)*sz
        data_R3[...,1] = sin(a)*sz
        data_R3[...,2] = cos(z)
        return data_R3
    @staticmethod
    def R3toS2(data_R3):
        """Transform R3 points to S2 azimuth, zenith angle pairs
            azimuth, zenith = arctan2(y,x), arctan2(hypot(x,y),z)
        """
        data_R3 = atleast_2d(data_R3)
        if data_R3.shape[-1] != 3:
            raise ValueError("data_R3 must have shape (...,3)")
        x,y,z = data_R3[...,0],data_R3[...,1],data_R3[...,2]
        data_S2 = empty(data_R3.shape[:-1]+(2,))
        data_S2[...,0] = arctan2(y,x)
        data_S2[...,1] = arctan2(hypot(x,y),z)
        return data_S2
    @staticmethod
    def CtoA(d):
        """unit-sphere chord length to arc-distance"""
        return 2*arcsin(d/2)
    @staticmethod
    def AtoC(a):
        """arcdistance to unit-sphere chord length"""
        return 2*sin(a/2)
    @staticmethod
    def arcdist(p0,p1):
        """Arcdistance between azimuth,zenith pairs
        
        Parameters
        ----------
        p0, p1 : array_like, shape (..., 2)
            p[...,0] = azimuth angles
            p[...,1] = zenith angles

        Returns
        -------
        ad :  array_like, shape is broadcast(p0,p1).shape
            Arcdistances between corresponding pairs in p0,p1
        """
        #formula comes from translating points into cartesian coordinates
        #taking the dot product to get the cosine between the two vectors
        #then arccos to return to angle, and simplify everything assuming real inputs
        p0,p1 = np.array(p0), np.array(p1)
        a0,z0 = p0[...,0], p0[...,1]
        a1,z1 = p1[...,0], p1[...,1]
        return np.arccos(np.cos(z0)*np.cos(z1)+np.cos(a0-a1)*np.sin(z0)*np.sin(z1))

class KNNRegressor:
    def __init__(self,x,y,k=10,tree_cls=None,y_mean=None,y_std=None,lazy=False):
        """KNNRegressor(x,y,k=10,tree_cls=None,y_mean=None,y_std=None)
        x,y : points, values
        k : number of y values to average per evaluation. tree leaf size is twice this
        tree_cls : class to use for tree. default is KDTree
        y_mean : function for taking mean of y values. default is np.mean
        y_std : function for taking std of y values. defautl is np.std
        """
        self._x = np.asarray(x) if x is not None else x
        self._y = np.asarray(y) if y is not None else y
        self._k = k
        self._ys = y_std
        self._ym = y_mean
        self._lazy = lazy
        if tree_cls is not None:
            self._tree_cls = tree_cls
        else:
            self._tree_cls = KDTree
        self._tree = None
        if not lazy:
            self.build_tree()
        if y_mean is None:
            self._ym = np.mean
        if y_std is None:
            self._ys = np.std

    def build_tree(self):
        if self._tree is None:
            self._tree = self._tree_cls(self._x,2*self._k)

    def query(self,x):
        self.build_tree()
        d,i = self._tree.query(x,self._k)
        return self._y[i]

    def dist(self,x):
        self.build_tree()
        d,i = self._tree.query(x,self._k)
        return np.mean(d,axis=-1)

    def value(self,x):
        return self._ym(self.query(x),axis=-1)

    def std(self,x):
        return self._ys(self.query(x),axis=-1)

    def stats(self,x):
        y = self.query(x)
        return self._ym(y,axis=-1), self._ys(y,axis=-1)

    def __call__(self,x):
        return self._ym(self.query(x),axis=-1)

    def merge(self,other):
        if self._x is None:
            nx = copy(other._x)
        elif other._x is None:
            return
        else:
            nx = concatenate((self._x,other._x),axis=0)
        if self._y is None:
            ny = copy(other._y)
        elif other._y is None:
            raise RuntimeError('Merged with bad KNNRegressor')
        else:
            ny = concatenate((self._y,other._y),axis=0)
        self._x = nx
        self._y = ny
        if not self._lazy:
            self.build_tree()

def S2KNNRegressor(x,y,k=10,lazy=False):
    return KNNRegressor(x,y,k,S2KDTree,angle_mean,angle_std,lazy)

def unused_fname(fname,sep='-',start=2,before_ext=True):
    if before_ext: root,ext = os.path.splitext(fname)
    else: root,ext = fname,''

    i,fn = start, fname
    while os.path.exists(fn):
        fn = root + sep + str(i) + ext
        i += 1
    return fn

# #invert the diverging colormaps:
# def invert_cmap_data(cmname,x=0.9):
#     return dict((c,[(v,1-x*a,1-x*b) for v,a,b in cm.datad[cmname][c]]) for c in ('red','green','blue'))

#Normalization function
def rescale(x,vmean=0.0,vstd=1.0,xmean=None,xstd=None):
    if xmean is None: xmean = mean(x)
    if xstd is None: xstd = std(x)
    return (x-xmean)*vstd/xstd + vmean

#def rescale(x,vmin=0.0,vmax=1.0,vrange=None,xmin=None,xmax=None):
#    if xmin is None: xmin = amin(x)
#    if xmax is None: xmax = amax(x)
#    return (x-xmin)*(vmax-vmin)/(xmax-xmin) + vmin

#Wrap???
def wrap(x,vmin=0.0,vmax=1.0):
    return (x-vmin)%(vmax-vmin)+vmin

#For calculating the Intensity, Dop, Aop.
def poincare(s):
    """return intensity, degree of polarization, angle of polarization"""
    #not actually the poincare sphere parameters, but whatever...
    s0 = array(s[...,0],copy=True)
    s1 = s[...,1]
    s2 = s[...,2]
    #hypot euclidean norm
    p = hypot(s1,s2)/s0 #degree of polarization
    a = arctan2(s2,s1)*0.5
    return (s0, p, a)

def az_zen_dist(p0,p1):
    """angular distance between two points in azimuth,zenith format"""
    #formula comes from translating points into cartesian coordinates
    #taking the dot product to get the cosine between the two vectors
    #then arccos to return to angle, and simplify everything assuming real inputs
    p0,p1 = array(p0), array(p1)
    a0,z0 = p0[...,0], p0[...,1]
    a1,z1 = p1[...,0], p1[...,1]
    return arccos(cos(z0)*cos(z1)+cos(a0-a1)*sin(z0)*sin(z1))

# Unused
def great_circle(lat0,lon0,head_deg,dist_deg):
    """return lat,lon dist_deg away from lat0,lon0"""
    ln = WGS84.Line(lat0,lon0,head_deg)
    #dist_deg - spherical arc length from the first point to the second in degrees
    #p - geodesic dictionary with keys lat1, lon1, azi1, lat2, lon2, azi2, s12, a12
    p = ln.ArcPosition(dist_deg)
    return p['lat2'],p['lon2']
# Unused
def geoline(lat0,lon0,az_deg,arcd_deg):
    """return array of lat,lon pairs starting from lat0,lon0 heading in az_deg, at arcd_deg away"""
    ln = WGS84.Line(lat0,lon0,az_deg)
    if len(shape(arcd_deg)) == 0:
        p = ln.ArcPosition(arcd_deg)
        return array((p['lat2'],p['lon2']))
    else:
        ad = array(arcd_deg)
        #Add one more dimension with degree of 2
        out = empty(shape(ad)+(2,))
        for i,a in enumerate(ad.flat):
            p = ln.ArcPosition(a)
            out[...,0].flat[i] = p['lat2']
            out[...,1].flat[i] = p['lon2']
        return out
# Unused
def load_file(fname,roi=()):
    with h5py.File(fname,'r') as f:
        roll = array(f['roll'])
        pitch = array(f['pitch'])
        heading = array(f['heading'])
        times = array(f['timestamps'])
        exps = array(f['exposures'])
        s = array(f['stokes'][roi])
        
        return roll,pitch,heading,s,times,exps

 # For calculating errors in two forms.
def sun_pos_error_l2(sun_head,sun_zen,cam_aop,cam_head,cam_pitch,m2):
    sim_aop = oceanaop(sun_head,sun_zen,cam_head,cam_pitch,m2)
    d = angle_diff(cam_aop,sim_aop,pi)
    return sum(d**2)

def sun_pos_error_l1(sun_head,sun_zen,cam_aop,cam_head,cam_pitch,m2):
    sim_aop = oceanaop(sun_head,sun_zen,cam_head,cam_pitch,m2)
    d = angle_diff(cam_aop,sim_aop,pi)
    return sum(abs(d))

#Rename function
#TODO:
sun_pos_error = sun_pos_error_l1
residuals_name = 'residuals-L1.h5'

def declination(lat,lon,elev,t,gm,radians=True):
    #TODO: CALCULATING THE DIFFERENCE BETWEEN TWO N
    b = broadcast(lat, lon, elev, t)
    res = empty(b.shape)
    for i,x in enumerate(b):
        la,lo,el,tt = x
        res.flat[i] = gm.declination(la,lo,el,gm.decimal_year(tt))
    if radians:
        res = deg2rad(res)
    return res

def sunpos_mag(t,lat,lon,elev,gm,temp=None,press=None,radians=True):
    #TODO: HOW TO COMPUTE SUNPOS HERE?
    """Observed sun heading, zenith -- using magnetic N"""
    #az_zen is a (...,5) dimension ndarray
    az_zen = sunpos(t,lat,lon,elev,temp,press,radians=radians)
    decl = declination(lat,lon,elev,t,gm,radians)
    az_zen[...,0] -= decl
    #subtract declination to go from true N to magnetic N
    return az_zen

def gps_error(lat,lon,t,sun_h,sun_z,gm):
    #can't pass sun_pos together because it doesn't vectorize right
    sun_pos = sun_h,sun_z
    fit_pos = sunpos_mag(t,lat,lon,0,gm)
    #use arc-distance
    return az_zen_dist(sun_pos,fit_pos)
gps_error_vec = vectorize(gps_error)


#No use!
def eval_model(t,lat,lon,head,pitch,tide=0,temp=None,press=None):
    """Model the sun position and AoP; 
        given:
          time, lat, lon,
          heading (radians, 0 is N, positive towards E),
          pitch (radians, 0 is horizontal),
          tide height (m above mean sea level),
          temperature (C), atmospheric pressure (mbar))
        returns:
          sun_az, sun_zen, sun_head (= az - decl), aop
    """
    #get the sun positions for each timestamp, at our known lat,lon
    #sun_head, sun_zen = sunpos_mag(t,lat,lon,tide,temp,press,radians=True)
    sun_head = sunpos_mag(t, lat, lon, tide, temp, press, radians=True)
    sun_zen = sun_head[...,1]
    sun_head = sun_head[...,0]

    #TODO: input and output argument mismatch
    #get the ocean model aop values for each camera position
    aop = oceanaop(sun_head,sun_zen,cam_head,cam_pitch,1.33)
    return sun_az,sun_zen,sun_head,aop

def fit_sun_head_to_aop(aop,head,pitch,ridx=1.33,verbose=False):
    #fit the sun heading to the data by assuming the sun zenith angle is pi/4:
    #TODO: the dimension of x0 should be the number of variables. So how many variables here?
    fit = minimize(sun_pos_error,x0=np.asarray([0]),args=(pi/4,aop,head,pitch,ridx),bounds=[(-2*pi,2*pi)],options={'gtol':1e-6})
    return fit.x

def fit_sun_to_aop(aop,head,pitch,sun_head_guess=None,ridx=1.33,verbose=False,vl=0):
    """fit_sun_to_aop, find sun position corresponding to aops at different headings & pitches, no time passing
    aop: array of AoPs, radians
    head: array of headings, radians, same shape as aop
    pitch: array of pitchs, radians, same shape as head
    ridx: refractive index of water
    """
    v = verbose
    _p(v,vl,'Fitting sun ',end='')
    if sun_head_guess is None:
        _p(v,0,'heading to AoP... ',end='',flush=True)
        #fit just the heading first to get a good initial point:
        sun_head_guess = fit_sun_head_to_aop(aop,head,pitch,ridx)
    _p(v,0,'heading and zenith to AoP... ',end='',flush=True)
    #now do both heading & zenith
    minfun = lambda x,*args: sun_pos_error(x[0],x[1],*args)
    #Originally x0=(sun_head_guess,pi/4)
    #TODO: verify modification
    fit = minimize(minfun,x0=np.asarray([sun_head_guess,pi/4]),args=(aop,head,pitch,ridx),bounds=[(-2*pi,2*pi),(.01,pi/2)],options={'gtol':1e-6})
    sun_hz = fit.x
    _p(v,0,'DONE',flush=True)
    return sun_hz

def fit_gps_to_sun(t,sun_pos,gm,lat_range=(-70,70),verbose=False,vl=0):
    v = verbose
    _p(v,vl,'Fitting gps to sun position... ',end='',flush=True)
    #start with grid search, 5 degree grid:
    lat_min,lat_max = min(lat_range),max(lat_range)
    #add: int() to the next two lines
    err_lat = linspace(lat_min,lat_max,int((lat_max-lat_min)/5 + 1))
    err_lon = linspace(0,355,int(355/5+1))
    err_lon,err_lat = meshgrid(err_lon,err_lat)
    #minimum error over grid:
    sun_h,sun_z = sun_pos
    min_idx = argmin(gps_error_vec(err_lat,err_lon,t,sun_h,sun_z,gm))
    x0 = (err_lat.flat[min_idx], err_lon.flat[min_idx])
    #search around that point now:
    minfun = lambda x, *args: gps_error(x[0],x[1],*args) #also log(gps_error(...))
    fit = minimize(minfun,x0=np.asarray([x0[0],x0[1]]),args=(t,sun_h,sun_z,gm),bounds=[(lat_min-5,lat_max+5),(-360,360)])
    fit_gps = fit.x
    _p(v,0,'DONE',flush=True)
    return fit_gps

def derivative(x,t):
    """take discrete derivative dx/dt and resample t at midpoints. returns (dx/dt, t_resampled). Operates on last axis"""
    dxdt = diff(x)/diff(t)
    t_r = (t[...,:-1] + t[...,1:])/2
    return dxdt,t_r

def sinc_filter(n,o):
    x = sinc(linspace(-o,o,n+2)[1:-1])
    return x/sum(x)

#No use!
def correct_imu(r,p,h,t,dt=0,n=8,decorr=False):
    """ correct_imu(r,p,h,t,dt=0,n=8,decorr=False)
    r,p,h,t: roll, pitch, heading, time
    dt: difference between t and when r,p,h was actually sampled
        sample time = t - dt
    n: number of filter taps = 2*n+1
    decorr: if true, remove components of r and p that correlate with 1st and 2nd derivatives of h

    Returns: r,p,h, dh
    """
    
    h = unwrap(h)
    #derivative:
    dh,tdh = derivative(h,t)
    #2nd derivative:
    if decorr: d2h,td2h = derivative(dh,tdh)
    #interpolate to sync to t
    #derivatives lost up to 2 points on either end, interp will just extend
    dh = interp(t,tdh,dh)
    if decorr: d2h = interp(t,td2h,d2h)
    #filter everything
    f = sinc_filter(2*n+1,3)
    tf = convolve(t,f,'valid')-dt #adjust time
    rf = convolve(r,f,'valid')
    pf = convolve(p,f,'valid')
    hf = convolve(h,f,'valid')
    dhf = convolve(dh,f,'valid')
    if decorr: d2hf = convolve(d2h,f,'valid')
    
    #remove components of roll, pitch that correlate with changes in heading
    # and interpolate filtered signals back to original sampling times
    if decorr:
        #least squares
        a = stack((dhf,d2hf,ones_like(d2hf)),-1)
        xr = lstsq(a,rf)[0]
        xp = lstsq(a,pf)[0]
        rc = interp(t,tf,rf - xr[0]*dhf - xr[1]*d2hf)
        pc = interp(t,tf,pf - xp[0]*dhf - xp[1]*d2hf)
    else:
        rc = interp(t,tf,rf)
        pc = interp(t,tf,pf)
    hc = interp(t,tf,hf)
    dhc = interp(t,tf,dhf)
    return rc,pc,hc,dhc

def file_replace(path0,path1):
    #replace path0 with path1
    if os.path.exists(path0):
        os.rename(path0,path0+'.old')
        os.rename(path1,path0)
        os.remove(path0+'.old')
        return True
    else:
        os.rename(path1,path0)
        return False

class KernelRegression:
    def __init__(self,sigma_deg=1,grid_size=2):
        self.sigma_deg = sigma_deg
        self.grid_size = grid_size

        self.sigma = deg2rad(self.sigma_deg)
        self.head = linspace(-pi,pi,ceil(2*pi/(self.grid_size*self.sigma)),endpoint=False)
        self._h = self.head % (2*pi)
        if self.sigma < 0.1:
            self.kernel = stats.norm(scale=self.sigma)
        else:
            self.kernel = stats.vonmises(self.sigma**-2)
        self.weights = None
        self.s1_weighted = None
        self.s2_weighted = None
        
        self._valid = False
        self._s1 = None
        self._s2 = None

    def _invalidate(self):
        self._valid = False
        self._s1=None
        self._s2=None

    def _apply_weights(self):
        if not self._valid:
            s1 = self.s1_weighted / self.weights
            s2 = self.s2_weighted / self.weights
            self._s1 = s1
            self._s2 = s2
            self._valid = True

    def do_regression(self,heads,aops,verbose=False,vl=0):
        v = verbose
        _p(v,vl,'Doing kernel regression... ',end='',flush=True)
        #evaluate kernel
        k = self.kernel.pdf(self.head[:,None] - heads[None,:]) # (h x head_rel)
        #convert angle residuals to stokes space so that we can take the mean
        s1,s2 = cos(2*aops), sin(2*aops)
        #use the kernel to create the weighted sums and total weights
        w = sum(k,axis=-1)
        s1w = sum(k*s1[None,:],axis=-1) #need to divide by w to get actual prediction
        s2w = sum(k*s2[None,:],axis=-1)
        self._invalidate()
        self.weights = w
        self.s1_weighted = s1w
        self.s2_weighted = s2w
        _p(v,0,'DONE',flush=True)

    def merge(self,other):
        if self.sigma_deg != other.sigma_deg or self.grid_size != other.grid_size:
            raise ValueError('KernelRegressions can\'t be merged: sigma_deg or grid_size do not match')
        if self.weights is not None:
            w = self.weights + other.weights
            s1w = self.s1_weighted + other.s1_weighted
            s2w = self.s2_weighted + other.s2_weighted
        else:
            w = other.weights
            s1w = other.s1_weighted
            s2w = other.s2_weighted
        self._invalidate()
        self.weights = w
        self.s1_weighted = s1w
        self.s2_weighted = s2w

    def interpolate(self,h):
        self._apply_weights()
        s1 = interp(h,self.head,self._s1,period=2*pi)
        s2 = interp(h,self.head,self._s2,period=2*pi)
        return 0.5*arctan2(s2,s1)

    def evaluate_spin(self,sv):
        """evaluate given a SpinVideo"""
        return self.interpolate(sv.head_rel)

    __call__ = evaluate_spin

class SpinVideo:
    def __init__(self,dname,calfile,roi,bad_frames,depth,gps,site,gm):
        self.gm = gm #geomagnetic model
        #video properties
        self.dir_name = dname #directory name
        self.name = os.path.basename(dname) #video name
        self.trace_path = os.path.join(dname,'trace.h5')
        self.model_path = os.path.join(dname,'model.h5')
        self.res_path = os.path.join(dname,residuals_name)
        self.calfile = calfile #calibration file
        self.roi = roi # r0, c0, h, w; r0,c0 are center of roi
        r0,c0,h,w = roi
        self.roi_slice = (slice(r0-h//2,r0+(h+1)//2),slice(c0-w//2,c0+(w+1)//2))
        self.bad_frames = bad_frames #list of ranges of bad frames
        self._bad_frames = bad_frames
        self.depth = depth
        self.gps = gps #lat,lon,elevation of site where data was taken
        self.site = site #name of site where data was taken

        #measurements
        self.time = None #timestamp; in seconds
        self.time_mean = None
        self.stokes = None #mean stokes parameters in roi
        self.roll = None #roll; in radians
        self.pitch = None #pitch; in radians
        self.head = None #magnetic heading; in radians
        #computed from measurements:
        self.intensity = None #total intensity (S0)
        self.angle = None #e-vector angle (AoP)
        self.partial = None #partial polarization (DoLP)
        self.head_speed = None #magnetic heading rate of change; in radians/second

        #model--based on t, head, pitch, gps
        self.sun_pos = None #actual sun magnetic heading, zenith; in radians
        self.sun_pos_mean = None #sun position at the mean timestamp. magnetic heading, zenith; in radians
        self.model_stokes = None #model stokes vector evaluated for head,pitch,sun_pos
        self.model_intensity = None
        self.model_angle = None
        self.model_partial = None
        #residuals--based on measurements and model
        self.sun_head_est = None #estimated sun magnetic heading; in radians
        self.head_rel = None #heading relative to sun: head - sun_head_est
        self.model_aop_res = None #cam_aop - model_aop
        #kernel regression of residuals
        self.regression = None
        self.model_aop_res_est = None #estimate of model_aop_res
        #position estimates
        self.fit_sun_pos = None
        self.fit_gps = None
        self.fit_sun_error = None #arcdistance
        self.fit_gps_error = None #meters

    def load_all(self,force=False,verbose=False,vl=0):
        self.get_measurements(force,verbose,vl)
        self.get_model(force,verbose,vl)
        self.get_residuals(force,verbose,vl)

    def _smooth_compass(self,n=8,dt=0):
        #unwrap heading to avoid large jumps
        h = unwrap(self.head)
        #take derivative to get speed
        dh,tdh = derivative(h,self.time)
        #interpolate to original times
        dh = interp(self.time,tdh,dh)
        #filter
        f = sinc_filter(2*n+1,3)
        tf = convolve(self.time,f,'valid')-dt #adjust time
        rf = convolve(self.roll,f,'valid')
        pf = convolve(self.pitch,f,'valid')
        hf = convolve(h,f,'valid')
        dhf = convolve(dh,f,'valid')
        #interpolate back to original times
        rc = interp(self.time,tf,rf)
        pc = interp(self.time,tf,pf)
        hc = interp(self.time,tf,hf)
        dhc = interp(self.time,tf,dhf)
        #store
        self.roll = rc
        self.pitch = pc
        self.head = hc
        self.head_speed = dhc

    def _clip_bad_frames(self):
        if self.bad_frames:
            #we need to mask out the bad frames
            #make a boolean mask of which frames are good
            mask = ones(self.time.shape[0],bool8) #all true to start
            for start,end in self.bad_frames:
                mask[start:end] = False #mark bad frames as false
            #find the indices of the good frames
            mask = find(mask)
            #index only the good frames out of the data arrays
            self.roll = self.roll[mask,...]
            self.pitch = self.pitch[mask,...]
            self.head = self.head[mask,...]
            self.head_speed = self.head_speed[mask,...]
            self.stokes = self.stokes[mask,...]
            self.intensity = self.intensity[mask,...]
            self.angle = self.angle[mask,...]
            self.partial = self.partial[mask,...]
            self.time = self.time[mask,...]
            self.bad_frames = []

    def get_measurements(self,force=False,verbose=False,vl=0):
        #load time, stokes, roll, pitch, heading from trace file
        #process the trace file if necessary
        v = verbose
        #don't load unless we need to:
        loaded = self.time is not None and self.stokes is not None
        loaded &= self.roll is not None and self.pitch is not None
        loaded &= self.head is not None
        if not force and loaded:
            return False
        _p(v,vl,'Loading {}... '.format(os.path.basename(self.trace_path)),end='',flush=True)
        if not os.path.exists(self.trace_path):
            _p(v,0,'DOES NOT EXIST',flush=True)
            self.process_measurements(v,vl+1)
            _p(v,vl,'Loading {}... '.format(os.path.basename(self.trace_path)),end='',flush=True)
        #check valid file:
        with h5py.File(self.trace_path,'r') as f:
            roi = array(f['roi'])
            bad_file = any(roi != array(self.roi))
        if bad_file:
            _p(v,0,'INVALID',flush=True)
            self.process_measurements(v,vl+1)
            _p(v,vl,'Loading {}... '.format(os.path.basename(self.trace_path)),end='',flush=True)
        with h5py.File(self.trace_path,'r') as f:
            t = array(f['timestamps'])
            s = array(f['stokes'])
            r = deg2rad(array(f['compass/roll']))
            p = deg2rad(array(f['compass/pitch']))
            h = deg2rad(array(f['compass/heading']))
            s0, dolp, aop = poincare(s)
            aop -= r
        _p(v,0,'DONE',flush=True)
        
        self.time = t
        self.time_mean = mean(t)
        self.stokes = s
        self.roll = r
        self.pitch = p
        self.head = h
        self.intensity = s0
        self.angle = aop
        self.partial = dolp
        self._smooth_compass() #also compute self.head_speed
        self._clip_bad_frames()
        return True

    def process_measurements(self,verbose=False,vl=0):
        v = verbose
        _p(v,vl,'Processing raw data of {}: {}...'.format(self.site,self.name),flush=True)
        vl += 1
        _p(v,vl,'Loading calibration file... ',end='',flush=True)

        gains,darks,pattern,slices = load_cal(self.calfile)
        _p(v,0,'DONE',flush=True)
        #only files that are named with 3 digits
        file_paths = sorted(glob(os.path.join(self.dir_name,'[0-9][0-9][0-9].h5')))
        #first we need to get the total length of the data
        nframes = 0
        cuts = [] #index where each file starts
        for fp in file_paths:
            with h5py.File(fp,'r') as f:
                cuts.append(nframes)
                nframes += len(f['timestamps'])
        #create a temporary output file
        tmp_path = self.trace_path+'.tmp'
        _p(v,vl,'Creating output file: {}'.format(os.path.basename(tmp_path)),flush=True)
        with h5py.File(tmp_path,'w') as outfile:
            #create datasets for processing metadata:
            outfile.create_dataset('cuts',data=array(cuts),dtype=np.int64)
            outfile.create_dataset('bad_frames',data=array(self.bad_frames),dtype=np.int64)
            outfile.create_dataset('roi',data=array(self.roi),dtype=np.int32)
            outfile.create_dataset('depth',data=array(self.depth),dtype=np.float32)
            outfile.create_dataset('gps',data=array(self.gps),dtype=np.float32)
            #bug of pycharm
            outfile.create_dataset('site',data=u'{}'.format(self.site),dtype=h5py.special_dtype(vlen=unicode))
            #create datasets for measurements:
            outfile.create_dataset('compass/accelX',shape=(nframes,),dtype=np.float32)
            outfile.create_dataset('compass/accelY',shape=(nframes,),dtype=np.float32)
            outfile.create_dataset('compass/accelZ',shape=(nframes,),dtype=np.float32)
            outfile.create_dataset('compass/calstatus',shape=(nframes,),dtype=np.int8)
            outfile.create_dataset('compass/distortion',shape=(nframes,),dtype=np.int8)
            outfile.create_dataset('compass/fresh',shape=(nframes,),dtype=np.int8)
            outfile.create_dataset('compass/heading',shape=(nframes,),dtype=np.float32)
            outfile.create_dataset('compass/magX',shape=(nframes,),dtype=np.float32)
            outfile.create_dataset('compass/magY',shape=(nframes,),dtype=np.float32)
            outfile.create_dataset('compass/magZ',shape=(nframes,),dtype=np.float32)
            outfile.create_dataset('compass/pitch',shape=(nframes,),dtype=np.float32)
            outfile.create_dataset('compass/roll',shape=(nframes,),dtype=np.float32)
            outfile.create_dataset('compass/temperature',shape=(nframes,),dtype=np.float32)
            outfile.create_dataset('exposures',shape=(nframes,),dtype=np.int64)
            outfile.create_dataset('timestamps',shape=(nframes,),dtype=np.float64)
            outfile.create_dataset('zooms',shape=(nframes,),dtype=np.int16)
            outfile.create_dataset('stokes',shape=(nframes,3),dtype=np.float64)

            #process each file in turn:
            for c,fp in zip(cuts,file_paths):
                _p(v,vl,'{}: '.format(os.path.basename(fp)),end='',flush=True)
                with h5py.File(fp,'r') as f:
                    n = len(f['timestamps'])

                    _p(v,0,'metadata... ',end='',flush=True)
                    exps = array(f['exposures'])
                    outfile['exposures'][c:c+n] = exps
                    outfile['timestamps'][c:c+n] = array(f['timestamps'])/1000 #convert to seconds
                    for ds in ['compass/accelX','compass/accelY','compass/accelZ','compass/calstatus',
                                'compass/distortion','compass/fresh','compass/heading','compass/magX',
                                'compass/magY','compass/magZ','compass/pitch','compass/roll',
                                'compass/temperature','zooms']:
                        outfile[ds][c:c+n] = f[ds][:]
                    
                    _p(v,0,'load... ',end='',flush=True)
                    data = squeeze(array(f['data']))
                    _p(v,0,'process... ',end='',flush=True)
                    #allocate stokes, process by blocks of frames
                    stokes = empty((n,3))
                    prog = 0
                    for i in range(0,data.shape[0],10):
                        while i > n*prog/10:
                            _p(v,0,'{}'.format(prog),end='',flush=True)
                            prog += 1
                        j = min(i+10,data.shape[0])
                        s = process(data[i:j],(gains,darks,pattern,slices),'none')
                        #take mean over ROI and normalize to exposure
                        s = s[(slice(None),)+self.roi_slice].reshape((s.shape[0],-1,s.shape[-1]))
                        stokes[i:j] = mean(s,axis=1) / exps[i:j,None]
                    del data
                    _p(v,0,' write... ',end='',flush=True)
                    outfile['stokes'][c:c+n] = stokes
                _p(v,0,'DONE',flush=True)
                #done with input file
            #done with loop over input files
        #done with temporary file
        _p(v,vl,'Renaming {} to {}'.format(os.path.basename(tmp_path),os.path.basename(self.trace_path)),flush=True)
        if file_replace(self.trace_path,tmp_path):
            _p(v,vl,'Replaced {}'.format(os.path.basename(self.trace_path)),flush=True)
        _p(v,vl-1,'+DONE',flush=True)

    def get_model(self,force=False,verbose=False,vl=0):
        v = verbose
        loaded = self.sun_pos is not None and self.sun_pos_mean is not None
        loaded &= self.model_stokes is not None
        if not force and loaded:
            return False
        
        _p(v,vl,'Loading {}... '.format(os.path.basename(self.model_path)),end='',flush=True)
        if not os.path.exists(self.model_path):
            _p(v,0,'DOES NOT EXIST',flush=True)
            self.process_model(v,vl+1)
            _p(v,vl,'Loading {}... '.format(os.path.basename(self.model_path)),end='',flush=True)
        if os.path.exists(self.trace_path) and file_newer(self.trace_path, self.model_path):
            _p(v,0,'INVALID',flush=True)
            self.process_model(v,vl+1)
            _p(v,vl,'Loading {}... '.format(os.path.basename(self.model_path)),end='',flush=True)
        with h5py.File(self.model_path,'r') as f:
            sp = array(f['sun_pos'])
            spm = array(f['sun_pos_mean'])
            ms = array(f['model_stokes'])
            s0, dolp, aop = poincare(ms)
        _p(v,0,'DONE',flush=True)
        self.sun_pos = sp
        self.sun_pos_mean = spm
        self.model_stokes = ms
        self.model_intensity = s0
        self.model_angle = aop
        self.model_partial = dolp
        return True

    def process_model(self,verbose=False,vl=0):
        v = verbose
        _p(v,vl,'Evaluating models for {}: {}...'.format(self.site,self.name),flush=True)
        vl += 1
        self.get_measurements(False,v,vl)
        _p(v,vl,'Evaluating sun position... ',end='',flush=True)
        lat,lon,elev = self.gps
        sun_pos = sunpos_mag(self.time,lat,lon,elev,self.gm,radians=True)
        sun_pos_mean = sunpos_mag(self.time_mean,lat,lon,elev,self.gm,radians=True)
        _p(v,0,'DONE',flush=True)
        _p(v,vl,'Evaluating ocean model... ',end='',flush=True)
        model_stokes = oceanstokes(sun_pos[...,0],sun_pos[...,1],self.head,self.pitch,1.33)
        _p(v,0,'DONE',flush=True)

        tmp_path = self.model_path+'.tmp'
        _p(v,vl,'Writing to {}...'.format(os.path.basename(tmp_path)),end='',flush=True)
        with h5py.File(tmp_path,'w') as f:
            f.create_dataset('sun_pos',data=sun_pos)
            f.create_dataset('sun_pos_mean',data=sun_pos_mean)
            f.create_dataset('model_stokes',data=model_stokes)
        _p(v,0,'DONE',flush=True)
        _p(v,vl,'Renaming {} to {}'.format(os.path.basename(tmp_path),os.path.basename(self.model_path)),flush=True)
        if file_replace(self.model_path,tmp_path):
            _p(v,vl,'Replaced {}'.format(os.path.basename(self.model_path)),flush=True)
        _p(v,vl-1,'+DONE',flush=True)

    def get_residuals(self,force=False,verbose=False,vl=0):
        v = verbose

        loaded = self.sun_head_est is not None and self.head_rel is not None
        loaded &= self.model_aop_res is not None
        if not force and loaded:
            return False

        _p(v,vl,'Loading {}... '.format(os.path.basename(self.res_path)),end='',flush=True)
        if not os.path.exists(self.res_path):
            _p(v,0,'DOES NOT EXIST',flush=True)
            self.process_residuals(v,vl+1)
            _p(v,vl,'Loading {}... '.format(os.path.basename(self.res_path)),end='',flush=True)
        if os.path.exists(self.res_path) and file_newer(self.model_path,self.res_path):
            _p(v,0,'INVALID',flush=True)
            self.process_residuals(v,vl+1)
            _p(v,vl,'Loading {}... '.format(os.path.basename(self.res_path)),end='',flush=True)
        with h5py.File(self.res_path) as f:
            she = array(f['sun_head_est'])
            hr = array(f['cam_head_rel'])
            mar = array(f['model_aop_res'])
        _p(v,0,'DONE',flush=True)
        self.sun_head_est = she
        self.head_rel = hr
        self.model_aop_res = mar
        return True

    def process_residuals(self,verbose=False,vl=0):
        v = verbose
        _p(v,vl,'Evaluating model residuals for {}: {}...'.format(self.site,self.name))
        vl += 1
        self.get_measurements(False,v,vl)
        self.get_model(False,v,vl)
        _p(v,vl,'Estimating sun heading from measurements... ',end='',flush=True)
        sun_head_est = fit_sun_head_to_aop(self.angle,self.head,self.pitch)
        _p(v,0,'DONE',flush=True)
        _p(v,vl,'Computing angle residuals... ',end='',flush=True)
        cam_head_rel = angle_diff(self.head,sun_head_est) #0 is approximately towards sun
        model_aop_res = angle_diff(self.angle,self.model_angle,pi)
        _p(v,0,'DONE',flush=True)
        tmp_path = self.res_path+'.tmp'
        _p(v,vl,'Writing to {}...'.format(os.path.basename(tmp_path)),end='',flush=True)
        with h5py.File(tmp_path) as f:
            f.create_dataset('sun_head_est',data=sun_head_est)
            f.create_dataset('cam_head_rel',data=cam_head_rel)
            f.create_dataset('model_aop_res',data=model_aop_res)
        _p(v,0,'DONE',flush=True)
        _p(v,vl,'Renaming {} to {}'.format(os.path.basename(tmp_path),os.path.basename(self.res_path)),flush=True)
        if file_replace(self.res_path,tmp_path):
            _p(v,vl,'Replaced {}'.format(os.path.basename(self.res_path)),flush=True)
        _p(v,vl-1,'+DONE',flush=True)

    def do_kernel_regression(self,sigma_deg=1,grid_size=2,verbose=False,vl=0):
        """do a kernel regression on model_aop_res vs head_rel, with grid spacing ~2*sigma_deg"""
        v = verbose
        _p(v,vl,'Kernel regression of angle residuals...',flush=True)
        self.get_residuals(False,v,vl+1)
        kr = KernelRegression(sigma_deg,grid_size)
        kr.do_regression(self.head_rel,self.model_aop_res)
        self.regression = kr
        _p(v,vl,'DONE')

    def estimate_sun(self,angle_res_func=None,verbose=False,vl=0):
        v = verbose
        _p(v,vl,'Estimating sun position...',flush=True)
        vl += 1
        self.get_measurements(False,v,vl)
        self.get_model(False,v,vl)
        self.get_residuals(False,v,vl) #not strictly necessary, but loads sun_head_est which is useful
        a_res_est = 0 #angle residual estimate
        if angle_res_func is not None:
            _p(v,vl,'Estimating model residuals... ',end='',flush=True)
            a_res_est = angle_res_func(self)
            self.model_aop_res_est = a_res_est
            _p(v,0,'DONE',flush=True)
        fit_sun_pos = fit_sun_to_aop(self.angle-a_res_est,self.head,self.pitch,self.sun_head_est,1.33,v,vl)
        fit_sun_error = az_zen_dist(self.sun_pos_mean,fit_sun_pos)
        self.fit_sun_pos = fit_sun_pos
        self.fit_sun_error = fit_sun_error
        _p(v,vl-1,'+DONE: {} error = {:.3f}\xB0'.format(fit_sun_pos,rad2deg(fit_sun_error)))

    def estimate_position(self,lat_range=(-70,70),verbose=False,vl=0):
        v = verbose
        _p(v,vl,'Estimating position...',flush=True)
        vl += 1
        fit_gps = fit_gps_to_sun(self.time_mean,self.fit_sun_pos,self.gm,lat_range,v,vl)
        lat,lon = self.gps[:2]
        fit_lat,fit_lon = fit_gps
        fit_gps_error = WGS84.Inverse(lat,lon,fit_lat,fit_lon)['s12']
        self.fit_gps = fit_gps
        self.fit_gps_error = fit_gps_error
        _p(v,vl-1,'+DONE: {} error = {:.3f} km'.format(fit_gps,fit_gps_error/1000))

#No use!
def fit_residuals_knn(video_list,roi_hw,gm,rdir):
    #make a function from (relative heading,pitch) -> (aop model residual)
    #heading is relative to the estimated sun heading given by fit_sun_head_to_aop
    cam_head_rels, cam_pitchs, aop_ress = [],[],[]
    for dname, calfile, rc0, cam_pos, bad_frames in video_list:
        r0,c0,roi = make_roi(rc0,roi_hw)
        cam_meas = load_measurements(dname,calfile,roi,bad_frames,verbose=True)
        t, cam_stokes, cam_exps, cam_roll, cam_pitch, cam_head, cam_head_speed = cam_meas
        cam_s0,cam_dolp,cam_aop = poincare(cam_stokes)

        sun_head_est, cam_head_rel, model_aop_res = load_residuals(dname,calfile,cam_pos,rc0,roi_hw,bad_frames,gm)

        cam_head_rels.append(cam_head_rel)
        cam_pitchs.append(cam_pitch)
        aop_ress.append(model_aop_res)
    cam_head_rels,cam_pitchs,aop_ress = map(concatenate,(cam_head_rels,cam_pitchs,aop_ress))
    cam_hz = stack((cam_head_rels,pi/2-cam_pitchs),-1)
    return KNNRegressor(cam_hz,aop_ress,10,S2KDTree,angle_mean,angle_std)

def parse_video_list(fnames,root_dir):
    rel_video_list=[]
    if isscalar(fnames):
        fnames = [fnames]
    for fname in fnames:
        with open(fname,encoding='utf-8') as f:
            rel_video_list.extend(yaml.safe_load_all(f))
    abs_video_list=[]
    for site,folder,calfile,rc0,videos in rel_video_list:
        for dname,depth,gps_pos,bad_frames in videos:
            dname = os.path.abspath(os.path.join(root_dir,folder,dname))
            abs_video_list.append((site,dname,calfile,rc0,depth,gps_pos,bad_frames))
    return abs_video_list

#No use!
def test_imu_delays():
    #parse video file list, videos are in Data folder
    print('Loading video lists...')
    video_list = parse_video_list(sys.argv[1:],'../../Data')
    
    if len(video_list) == 0:
        print('No videos to process!')
        sys.exit()
    
    #print which videos we'll be processing:
    for dname,calfile,cam_pos,bad_frames in video_list:
        print('{} ({}):\n    GP = {}, {}; D = {}; E = {}'.format(dname,calfile,*cam_pos))

    #set up geomagnetic model, but don't load data until necessary
    gm = EMMMesh('geomag/data/EMM-720_V3p1_static.bin','geomag/data/EMM-720_V3p1_secvar.bin', delay_load = True)

    #results directory
    rdir = unused_fname(datetime.now().strftime('run-%Y-%m-%d-%H.%M.%S'),before_ext=False)
    os.mkdir(rdir)
    print('Saving results to {}'.format(rdir))

    #geometric camera calibration
    r0,c0 = 495//2,776//2 #center
    roiw, roih = 100,100 #width & height
    roi = (slice(r0-roih//2,r0+(roih+1)//2),slice(c0-roiw//2,c0+(roiw+1)//2))
    print('Using ROI = ',roi)

    #process the data
    summary_data = []
    for dname, calfile, cam_pos, bad_frames in video_list:
        for imu_dt in linspace(0,1000,11):
            sd = estimate_position(dname,calfile,cam_pos,(r0,c0),roi,bad_frames,gm,rdir,imu_dt)
            if sd is not None:
                summary_data.append(sd)

    if len(summary_data) == 0:
        print('All failed!')
        sys.exit(1)
        
    print('Saving summary data... ',end='',flush=True)
    summary_data = list(map(array,zip(*summary_data))) #convert to arrays
    summary_data[0] = summary_data[0].astype('S') #fixed length strings for the video names
    ds_names = 'name','t_mean','cam_lat','cam_lon','cam_depth','cam_decl','sun_head_mean','sun_zen_mean','fit_head','fit_zen','fit_sun_error','fit_lat','fit_lon','fit_decl','fit_gps_error'
    with h5py.File(os.path.join(rdir,'summary.h5'),'w') as f:
        for n,d in zip(ds_names,summary_data):
            f.create_dataset(n,data=d)
    print('DONE',flush=True)

def plot_map(sv,sfx='',rdir=None):
    fig = figure(sv.name+sfx)
    lat,lon,elev = sv.gps
    rlat,rlon = round(lat,-1),round(lon,-1) #to nearest 10°
    fit_lat,fit_lon = sv.fit_gps
    ## Grid for error map
    err_lat = linspace(lat-20,lat+20,81)
    err_lon = linspace(lon-25,lon+25,101)
    err_lon,err_lat = meshgrid(err_lon,err_lat)
    #heading map, zenith map
    hz_map = sunpos_mag(sv.time_mean,err_lat,err_lon,0,sv.gm)
    head_map, zen_map = hz_map[...,0], hz_map[...,1]
    #heading errors, zenith errors
    head_err_map = angle_diff(head_map,sv.sun_pos_mean[0]) # fit_head
    zen_err_map = zen_map - sv.sun_pos_mean[1] # fit_zen
    el_err_map = -zen_err_map
    #arcdistance map
    #arc_err_map = gps_error_vec(err_lat,err_lon,t_mean,fit_sun_pos,gm)

    mp = Basemap(lon-25,lat-20,lon+25,lat+20,resolution='i')
    mp.drawmapboundary(fill_color='#bbf2f2')
    mp.fillcontinents(lake_color='#bbf2f2')

    mp.drawparallels(rlat+linspace(-30,30,7),dashes=[10,10],labels=[1,1,1,1],zorder=0,color='#99cccc')
    mp.drawmeridians(rlon+linspace(-30,30,7),dashes=[10,10],labels=[1,1,1,1],zorder=0,color='#99cccc')
    #mp.imshow(arc_err_map,cmap=cm.jet,interpolation='nearest')
    mp.contour(err_lon,err_lat,rad2deg(head_err_map),levels=linspace(-10,10,5),cmap='CyKOr',latlon=True)
    mp.contour(err_lon,err_lat,head_err_map,levels=[angle_diff(sv.fit_sun_pos[0],sv.sun_pos_mean[0])],colors='r',linestyles=':',latlon=True)
    mp.contour(err_lon,err_lat,rad2deg(el_err_map),levels=linspace(-10,10,5),cmap='CyKOr',latlon=True)
    mp.contour(err_lon,err_lat,zen_map,levels=[sv.fit_sun_pos[1]],colors='r',linestyles=':',latlon=True)
    #mp.drawgreatcircle(lon,lat,cam_to_sun[1],cam_to_sun[0],color='g')
    mp.plot(lon,lat,'gD',latlon=True)
    #mp.drawgreatcircle(fit_lon,fit_lat,fit_to_sun[1],fit_to_sun[0],color='r')
    mp.plot(fit_lon,fit_lat,'rs',latlon=True)

    xlabel('Longitude'); ylabel('Latitude')
    title(sv.name+'''
        True Lat, Lon = {la:.5f}\xB0, {lo:.5f}\xB0
        Fit Lat, Lon = {fla:.5f}\xB0, {flo:.5f}\xB0; Error = {er:.2f} mi'''.format(la=lat,lo=lon,fla=fit_lat,flo=fit_lon,er=units.convert(sv.fit_gps_error,'m','mi')))
    tight_layout()
    if rdir is not None:
        savefig(unused_fname(os.path.join(rdir,sv.name+sfx+'.png')),dpi=100)
        close(fig)
    return fig

def plot_map_multi(svs,name,rdir=None):
    #assumes that all svs have the same base location
    
    fig = figure(name+'-map')
    lat,lon,elev = svs[0].gps
    rlat,rlon = round(lat,-1),round(lon,-1) #to nearest 10°
    mp = Basemap(lon-25,lat-20,lon+25,lat+20,resolution='i')
    mp.drawmapboundary(fill_color='#bbf2f2')
    mp.fillcontinents(lake_color='#bbf2f2')

    mp.drawparallels(rlat+linspace(-30,30,7),dashes=[10,10],labels=[1,1,1,1],zorder=0,color='#99cccc')
    mp.drawmeridians(wrap(rlon+linspace(-30,30,7),-180,180),dashes=[10,10],labels=[1,1,1,1],zorder=0,color='#99cccc')
    #true location
    mp.plot(lon,lat,'gD',latlon=True)
    for sv in svs:
        fit_lat,fit_lon = sv.fit_gps
    
        mp.plot(fit_lon,fit_lat,'rs',latlon=True)

    if rdir is not None:
        savefig(unused_fname(os.path.join(rdir,name+'-map.png')),dpi=100)
        close(fig)
    return fig

def write_csv(svs,fname):
    fields = ['name','site','time','lat','lon','depth','sun_head','sun_zen',
            'sun_head_est','sun_zen_est','sun_pos_error','lat_est','lon_est','gps_error']
    with open(fname,'w') as f:
        f.write(', '.join(fields))
        f.write('\n')
        for sv in svs:
            name = sv.name
            site = sv.site
            time = sv.time_mean
            depth = sv.depth
            lat,lon,elev = sv.gps
            sun_head, sun_zen = sv.sun_pos_mean
            sun_head_est, sun_zen_est = sv.fit_sun_pos
            sun_pos_error = sv.fit_sun_error
            lat_est, lon_est = sv.fit_gps
            gps_error = sv.fit_gps_error
            f.write(', '.join(str(v) for v in [name,site,time,lat,lon,depth,sun_head,sun_zen,
                sun_head_est,sun_zen_est,sun_pos_error,lat_est,lon_est,gps_error]))
            f.write('\n')

def make_svs(video_list_file, roi_hw, gm):
    #svs = make_svs(video_list_file, [100, 100], gm)
    svs = []
    video_list = parse_video_list(video_list_file,'../../Data')
    # I don't know the itert procedure of this 'for' interation
    for site,dname,calfile,rc0,depth,gps_pos,bad_frames in video_list:
        roi = rc0 + roi_hw
        print('Loading {}'.format(dname),flush=True)
        svs.append(SpinVideo(dname,calfile,roi,bad_frames,depth,gps_pos,site,gm))
    return svs

def prep(video_list_file,svs,rdir,gm,mkdir=True):
    #prep(sys.argv[1:],None,None,gm,False)
    if svs is None and video_list_file is not None:
        print('Loading video lists...',end='',flush=True)
        svs = make_svs(video_list_file,[100,100],gm)

    if not svs:
        raise RuntimeError('No videos to process!')
    
    if rdir is None:
        rdir = unused_fname(datetime.now().strftime('../run-%Y-%m-%d-%H.%M.%S'),before_ext=False)
    
    if mkdir and not os.path.exists(rdir):
        os.mkdir(rdir)
    return svs,rdir

def preprocess(video_list_file=None,svs=None,gm=None,rdir=None):
    svs,rdir = prep(video_list_file,svs,rdir,gm)

    print('Saving results to {}'.format(rdir))
    summary_path = os.path.join(rdir,'traces.h5')
    residuals_path = os.path.join(rdir,'residuals.h5')

    print('Creating trace summary file {}'.format(os.path.basename(summary_path)))
    head_rel = []
    aop_res = []
    h5opts = {'shuffle':True,'compression':7}
    with h5py.File(summary_path,'w') as sumfile:
        #region of interest is 100x100:
        roi_hw = [100,100]
        #process each video:
        for sv in svs:
            print('Processing {}: {}'.format(sv.site, sv.name),flush=True)
            sv.load_all(verbose=True,vl=1)
            head_rel.append(sv.head_rel)
            aop_res.append(sv.model_aop_res)
            _p(1,1,'Adding data to {}... '.format(os.path.basename(summary_path)),end='',flush=True)
            grp = sumfile.create_group(sv.name)
            grp.create_dataset('name',data=sv.name) #scalars and small datasets don't need compression
            grp.create_dataset('calfile',data=sv.calfile)
            grp.create_dataset('roi',data=sv.roi)
            grp.create_dataset('depth',data=sv.depth)
            grp.create_dataset('gps',data=sv.gps)
            grp.create_dataset('site',data=sv.site)
            grp.create_dataset('time',data=sv.time,**h5opts)
            grp.create_dataset('stokes',data=sv.stokes,**h5opts)
            grp.create_dataset('roll',data=sv.roll,**h5opts)
            grp.create_dataset('pitch',data=sv.pitch,**h5opts)
            grp.create_dataset('head',data=sv.head,**h5opts)
            grp.create_dataset('intensity',data=sv.intensity,**h5opts)
            grp.create_dataset('angle',data=sv.angle,**h5opts)
            grp.create_dataset('partial',data=sv.partial,**h5opts)
            grp.create_dataset('head_speed',data=sv.head_speed,**h5opts)
            grp.create_dataset('sun_pos',data=sv.sun_pos,**h5opts)
            grp.create_dataset('sun_pos_mean',data=sv.sun_pos_mean)
            grp.create_dataset('model_stokes',data=sv.model_stokes,**h5opts)
            grp.create_dataset('model_intensity',data=sv.model_intensity,**h5opts)
            grp.create_dataset('model_angle',data=sv.model_angle,**h5opts)
            grp.create_dataset('model_partial',data=sv.model_partial,**h5opts)
            grp.create_dataset('sun_head_est',data=sv.sun_head_est)
            grp.create_dataset('head_rel',data=sv.head_rel,**h5opts)
            grp.create_dataset('model_aop_res',data=sv.model_aop_res,**h5opts)
            del grp
            _p(1,0,'DONE',flush=True)
    print('Creating residuals file {}... '.format(os.path.basename(residuals_path)),end='',flush=True)
    with h5py.File(residuals_path,'w') as resfile:
        head_rel = concatenate(head_rel)
        aop_res = concatenate(aop_res)
        resfile.create_dataset('head_rel',data=head_rel,**h5opts)
        resfile.create_dataset('model_aop_res',data=aop_res,**h5opts)
    print('DONE',flush=True)
    return svs

def reload_estimates(estimates_file,video_list_file=None,svs=None,gm=None):
    svs,rdir = prep(video_list_file,svs,None,gm,False)

    with h5py.File(estimates_file,'r') as estfile:
        for sv in svs:
            if sv.name in estfile:
                sv.fit_sun_pos = array(estfile[sv.name]['fit_sun_pos'])
                sv.fit_gps = array(estfile[sv.name]['fit_gps'])
                sv.fit_sun_error = array(estfile[sv.name]['fit_sun_error'])
                sv.fit_gps_error = array(estfile[sv.name]['fit_gps_error'])
    return svs

def estimate_no_residuals(video_list_file=None,svs=None,gm=None,rdir=None):
    svs,rdir = prep(video_list_file,svs,rdir,gm)

    print('Saving results to {}'.format(rdir))
    summary_path = os.path.join(rdir,'estimates-nr.h5')
    print('Creating trace summary file {}'.format(os.path.basename(summary_path)))
    with h5py.File(summary_path,'w') as sumfile:
        #process each video:
        for sv in svs:
            print('Processing {}: {}'.format(sv.site,sv.name),flush=True)
            # What is the function of syntax here?
            sv.estimate_sun(verbose=True,vl=1)
            sv.estimate_position(verbose=True,vl=1)
            _p(1,1,'Adding data to {}... '.format(os.path.basename(summary_path)),end='',flush=True)
            grp = sumfile.create_group(sv.name)
            grp.create_dataset('fit_sun_pos',data=sv.fit_sun_pos)
            grp.create_dataset('fit_gps',data=sv.fit_gps)
            grp.create_dataset('fit_sun_error',data=sv.fit_sun_error)
            grp.create_dataset('fit_gps_error',data=sv.fit_gps_error)
            del grp
            _p(1,0,'DONE',flush=True)
    return svs

# It is weired that this function has not been used.
def estimate_loo_kr(video_list_file=None,svs=None,sigma=1,grid=2,gm=None,rdir=None):
    svs,rdir = prep(video_list_file,svs,rdir,gm)

    print('Saving results to {}'.format(rdir))
    summary_path = os.path.join(rdir,'estimates_loo_kr.h5')
    print('Creating summary file {}'.format(os.path.basename(summary_path)))
    with h5py.File(summary_path,'w') as sumfile:
        sumfile.create_dataset('kr_params',data=array([sigma,grid],dtype=np.float64))
        #load the data
        krs = {} #kernel regressions, by site
        loo_krs = {} #leave-one-out kernel regressions, by site
        for sv in svs:
            print('Loading {}...'.format(sv.name),flush=True)
            if sv.site not in krs:
                krs[sv.site] = KernelRegression(sigma,grid)
                loo_krs[sv.site] = KernelRegression(sigma,grid)
            #load data
            sv.load_all(verbose=True,vl=1)
            #do Kernel regression on residuals
            sv.do_kernel_regression(sigma_deg=sigma,grid_size=grid,verbose=True,vl=1)
            krs[sv.site].merge(sv.regression)
        #make leave-one-out regressions
        for not_site in loo_krs:
            for site in krs:
                if site == not_site: continue
                loo_krs[not_site].merge(krs[site])
        #estimate positions using leave-one-out regressions
        for sv in svs:
            print('Processing {}: {}'.format(sv.site,sv.name),flush=True)
            sv.estimate_sun(verbose=True,vl=1)
            sv.estimate_position(loo_krs[sv.site],verbose=True,vl=1)
            _p(1,1,'Adding data to {}... '.format(os.path.basename(summary_path)),end='',flush=True)
            grp = sumfile.create_group(sv.name)
            grp.create_dataset('fit_sun_pos',data=sv.fit_sun_pos)
            grp.create_dataset('fit_gps',data=sv.fit_gps)
            grp.create_dataset('fit_sun_error',data=sv.fit_sun_error)
            grp.create_dataset('fit_gps_error',data=sv.fit_gps_error)
            del grp
            _p(1,0,'DONE',flush=True)
    return svs,krs,loo_krs

def estimate_loo_knn(video_list_file=None,svs=None,estimates_file=None,k=20,grid=1,gm=None,rdir=None):
    svs,rdir = prep(video_list_file,svs,rdir,gm)

    print('Saving results to {}'.format(rdir))
    summary_path = os.path.join(rdir,'estimates_loo_knn.h5')
    print('Creating summary file {}'.format(os.path.basename(summary_path)))
    with h5py.File(summary_path,'w') as sumfile:
        sumfile.create_dataset('knn_params',data=array([k,grid],dtype=np.float64))
        #load the data
        knns = {} #KNN regressions, by site
        loo_knns = {} #leave-one-out KNN regressions, by site
        for sv in svs:
            print('Loading {}...'.format(sv.name),flush=True)
            if sv.site not in knns:
                knns[sv.site] = S2KNNRegressor(None,None,k,True)
                loo_knns[sv.site] = S2KNNRegressor(None,None,k,True)
            #load data
            sv.load_all(verbose=True,vl=1)
        #load sun position estimates, compute others
        if estimates_file is not None:
            reload_estimates(estimates_file,svs=svs)
        #make knn regressors, compute any sun position estimates if necessary
        for sv in svs:
            if sv.fit_sun_pos is None:
                sv.estimate_sun(verbose=True,vl=1)
            h = angle_diff(sv.head, sv.fit_sun_pos[0]) #0 is approximately towards sun
            z = ones_like(h)*sv.fit_sun_pos[1]
            hz = stack((h,z),-1)
            knns[sv.site].merge(S2KNNRegressor(hz,sv.model_aop_res,k,True))

        #make leave-one-out knn regressors
        for not_site in loo_knns:
            for site in knns:
                if site == not_site: continue
                loo_knns[not_site].merge(knns[site])
        #sample the regressors on the grid:
        #h = linspace(-pi,pi,ceil(360/grid),endpoint=False)
        #z = linspace(deg2rad(5),pi/2,ceil(85/grid)) #leave out 5 degrees around the pole
        #hz = stack(meshgrid(h,z),-1)
        #loo_knns_sampled = {s:loo_knns[s].stats(hz) for s in loo_knns}
        def get_model_residuals(sv):
            #relative heading of video
            h = angle_diff(sv.head, sv.fit_sun_pos[0])
            z = ones_like(h)*sv.fit_sun_pos[1]
            hz = stack((h,z),-1)
            return loo_knns[sv.site].value(hz)

        #estimate positions using leave-one-out regressions
        for sv in svs:
            print('Processing {}: {}'.format(sv.site,sv.name),flush=True)
            _p(1,1,'Saving knn residuals...',end='',flush=True)
            with h5py.File(os.path.join(sv.dir_name,'residuals-knn-{}.h5'.format(k)),'w') as f:
                sv.model_aop_res_est = get_model_residuals(sv)
                f.create_dataset('model_aop_res_est',data=sv.model_aop_res_est)
            _p(1,0,'DONE',flush=True)
            sv.estimate_sun(get_model_residuals,verbose=True,vl=1)
            sv.estimate_position(verbose=True,vl=1)
            _p(1,1,'Adding data to {}... '.format(os.path.basename(summary_path)),end='',flush=True)
            grp = sumfile.create_group(sv.name)
            grp.create_dataset('fit_sun_pos',data=sv.fit_sun_pos)
            grp.create_dataset('fit_gps',data=sv.fit_gps)
            grp.create_dataset('fit_sun_error',data=sv.fit_sun_error)
            grp.create_dataset('fit_gps_error',data=sv.fit_gps_error)
            del grp
            _p(1,0,'DONE',flush=True)
    return svs,knns,loo_knns

def save_residuals(svs,estimates_file,fname):
    #write h5 file with residuals in a flat structure
    #include true sun position, estimated sun position, and relative camera heading for each
    
    #get number of data points first:
    n = 0
    for sv in svs:
        print('Loading {}... '.format(sv.name),flush=True)
        sv.load_all(verbose=True,vl=1)
        n += sv.model_aop_res.shape[0]
    print('Loading estimates... ',flush=True)
    reload_estimates(estimates_file,svs=svs)
    print('Saving residuals... ',flush=True)
    with h5py.File(fname,'w') as f:
        f.create_dataset('model_aop_res',dtype=float64,shape=(n,))
        f.create_dataset('cam_head_rel_est',dtype=float64,shape=(n,))
        f.create_dataset('cam_head_rel',dtype=float64,shape=(n,))
        f.create_dataset('sun_head',dtype=float64,shape=(n,))
        f.create_dataset('sun_zen',dtype=float64,shape=(n,))
        f.create_dataset('sun_head_est',dtype=float64,shape=(n,))
        f.create_dataset('sun_zen_est',dtype=float64,shape=(n,))
        i = 0
        for sv in svs:
            j = i+sv.model_aop_res.shape[0]
            f['sun_head'][i:j] = sv.sun_pos[...,0]
            f['sun_zen'][i:j] = sv.sun_pos[...,1]
            f['sun_head_est'][i:j] = sv.fit_sun_pos[0]
            f['sun_zen_est'][i:j] = sv.fit_sun_pos[1]
            f['cam_head_rel'][i:j] = sv.head - sv.sun_pos[...,0]
            f['cam_head_rel_est'][i:j] = sv.head - sv.fit_sun_pos[0]
            f['model_aop_res'][i:j] = sv.model_aop_res
            i = j

#No use!
def plot_residuals(res_file):
    #load data
    with h5py.File(res_file,'r') as f:
        model_aop_res = array(f['model_aop_res'])
        cam_head_rel = array(f['cam_head_rel'])
        cam_head_rel_est = array(f['cam_head_rel_est'])
        sun_head = array(f['sun_head'])
        sun_head_est = array(f['sun_head_est'])
        sun_zen = array(f['sun_zen'])
        sun_zen_est = array(f['sun_zen_est'])
    
    #plot vs true sun position
    norm = mpl.colors.Normalize(-30,30)
    residual_colors = apply_cmap(rad2deg(model_aop_res),cm.RdBu_r,norm)
    figure()
    for h,e,c in zip(wrap(rad2deg(cam_head_rel),-180,180),90-rad2deg(sun_zen),residual_colors):
        plot(h,e,'.',color=c,ms=3.0)
    xlim(-180,180); mxticks(labels='{}°',nticks=9)
    ylim(0,70); myticks(labels='{}°',nticks=8)
    xlabel('Heading relative to sun')
    ylabel('Sun elevation')
    make_colorbar(cm.RdBu_r,-30,30,ticks=linspace(-30,30,5),format='%1.0f°',label='AoP Model Error')
    tight_layout()

#No use!
def plot_trace_vs_model(svs):
    fig,axs = subplots(4,1)
    for ax,sv in zip(axs,svs):
        head = rad2deg(sv.head)
        aop = rad2deg(sv.angle)
        model_aop = rad2deg(sv.model_angle)
        sun_head = rad2deg(sv.sun_pos_mean[0])
        sun_elev = 90-rad2deg(sv.sun_pos_mean[1])
        mids,model_aop_bins = bin_wrap(head,model_aop,180,0,360,0,1,False)
        model_aop_mean = array(list(map(mean,model_aop_bins)))
        ax.plot(wrap(head,0,360),aop,'k.',ms=1)
        ax.plot(mids,model_aop_mean,'k')
        ax.plot([sun_head]*2,[-90,90],'k:')
        ax.set_xlim(0,360)
        ax.set_ylim(-70,70)
        ax.set_aspect('equal','box')
        mxticks(linspace(0,360,5),'',axes=ax)
        myticks(linspace(-45,45,3),'{:.0f}°',axes=ax)
    mxticks(linspace(0,360,5),'{:.0f}°',axes=ax)
    fig.tight_layout(h_pad=0)
    return fig,axs

#No use!
def main():
    #parse video file list, videos are in Data folder
    print('Loading video lists...')
    if len(sys.argv) > 2:
        training_list = parse_video_list(sys.argv[1],'../../Data')
        video_list = parse_video_list(sys.argv[2:],'../../Data')
    else:
        training_list = []
        video_list = parse_video_list(sys.argv[1],'../../Data')
    
    if len(video_list) == 0:
        print('No videos to process!')
        sys.exit()
    #print which videos we'll be processing:
    if len(training_list):
        print('VIDEOS FOR RESIDUAL ESTIMATION:')
        for dsite,dname,calfile,rc0,depth,gps_pos,bad_frames in training_list:
            lat,lon,e = gps_pos
            d = depth
            r0,c0 = rc0
            print('{}:\n GP = {}, {}; D = {}; E = {}; CTR = {}, {}; CAL = {}'.format(dname,lat,lon,d,e,r0,c0,calfile))
    print('VIDEOS TO PROCESS:')
    for dsite,dname,calfile,rc0,depth,gps_pos,bad_frames in video_list:
        lat,lon,e = gps_pos
        d = depth
        r0,c0 = rc0
        print('{}:\n GP = {}, {}; D = {}; E = {}; CTR = {}, {}; CAL = {}'.format(dname,lat,lon,d,e,r0,c0,calfile))

    #set up geomagnetic model, but don't load until necessary
    gm = EMMMesh('geomag/data/EMM-720_V3p1_static.bin','geomag/data/EMM-720_V3p1_secvar.bin',delay_load=True)

    #results directory
    rdir = unused_fname(datetime.now().strftime('run-%Y-%m-%d-%H.%M.%S'),before_ext=False)
    os.mkdir(rdir)
    print('Saving results to {}'.format(rdir))

    #estimate residuals:
    if len(training_list):
        print('Estimating residuals...',flush=True)
        aop_reg = fit_residuals_kr(training_list,(100,100),gm,rdir,deg2rad(1))
        print('DONE',flush=True)

        h = deg2rad(linspace(-180,180,720))
        aopr = aop_reg(h)
        fig=figure('AoP Model Residuals')
        plot(rad2deg(h),rad2deg(aopr))
        xlim(-180,180)
        mxticks(labels='{}\xB0'); myticks(labels='{}\xB0')
        tight_layout()
        savefig(os.path.join(rdir,'model-residual-pts.png'),dpi=100)
        close(fig)
        #make some figures:
        h,z = meshgrid(deg2rad(linspace(-180,180,720)),deg2rad(90-linspace(-5,25,120)))
        hz_grid = stack((h,z),-1)
        aopr = aop_reg(hz_grid)
        aops = aop_reg.std(hz_grid)

        fig=figure('AoP Model Residuals')
        for hz,ar in zip(rad2deg(aop_reg._x),rad2deg(aop_reg._y)):
            plot(hz[0],90-hz[1],'.',color=BuKRd((ar+15)/30))
        xlim(-180,180); ylim(-5,25)
        mxticks(labels='{}\xB0'); myticks(labels='{}\xB0')
        xlabel('Camera heading, relative to sun heading estimate')
        ylabel('Camera pitch')
        cbar=make_colorbar(BuKRd,-15,15)
        cbar.set_label('Residuals (deg)')
        tight_layout()
        savefig(os.path.join(rdir,'model-residual-pts.png'),dpi=100)
        close(fig)

        fig=figure('AoP Model Residuals')
        for hz,ar in zip(rad2deg(aop_reg._x),rad2deg(aop_reg._y)):
            plot(hz[0],90-hz[1],',',color=BuKRd((ar+15)/30))
        imshow(rad2deg(aopr),vmin=-15,vmax=15,cmap=BuKRd,aspect='auto',origin='lower',extent=(-180,180,-5,25))
        xlim(-180,180); ylim(-5,25)
        mxticks(labels='{}\xB0'); myticks(labels='{}\xB0')
        xlabel('Camera heading, relative to sun heading estimate')
        ylabel('Camera pitch')
        cbar=make_colorbar(BuKRd,-15,15)
        cbar.set_label('Estimated Residuals (deg)')
        tight_layout()
        savefig(os.path.join(rdir,'model-residual-reg.png'),dpi=100)
        close(fig)

        fig=figure('AoP Model Residuals')
        for hz,ar in zip(rad2deg(aop_reg._x),rad2deg(aop_reg._y)):
            plot(hz[0],90-hz[1],',',color=BuKRd((ar+15)/30))
        imshow(rad2deg(aops),vmin=0,vmax=5,cmap=cm.viridis,aspect='auto',origin='lower',extent=(-180,180,-5,25))
        xlim(-180,180); ylim(-5,25)
        mxticks(labels='{}\xB0'); myticks(labels='{}\xB0')
        xlabel('Camera heading, relative to sun heading estimate')
        ylabel('Camera pitch')
        colorbar().set_label('Std (deg)')
        cbar=make_colorbar(BuKRd,-15,15)
        cbar.set_label('Residuals (deg)')
        tight_layout()
        savefig(os.path.join(rdir,'model-residual-std.png'),dpi=100)
        close(fig)
        return
        del cam_hz,aop_res
    else:
        aop_res_func = None
    #process the data
    summary_data = []
    for dsite,dname,calfile,rc0,depth,gps_pos,bad_frames in video_list:

        sd = estimate_position(dname,calfile,cam_pos,rc0,(100,100),bad_frames,aop_res_func,gm,rdir)
        if sd is not None:
            summary_data.append(sd)

    if len(summary_data) == 0:
        print('All failed!')
        sys.exit(1)
        
    print('Saving summary data... ',end='',flush=True)
    summary_data = list(map(array,zip(*summary_data))) #convert to arrays
    summary_data[0] = summary_data[0].astype('S') #fixed length strings for the video names
    ds_names = 'name','t_mean','cam_lat','cam_lon','cam_depth','cam_decl','sun_head_mean','sun_zen_mean','fit_head','fit_zen','fit_sun_error','fit_lat','fit_lon','fit_decl','fit_gps_error'
    with h5py.File(os.path.join(rdir,'summary.h5'),'w') as f:
        for n,d in zip(ds_names,summary_data):
            f.create_dataset(n,data=d)
    print('DONE',flush=True)


# Main function
if __name__ == '__main__':
    try:
        #set up geomagnetic model, but don't load until necessary
        gm = EMMMesh('geomag/data/EMM-720_V3p1_static.bin','geomag/data/EMM-720_V3p1_secvar.bin',delay_load=True)

        svs,rdir = prep(sys.argv[1:],None,None,gm,False)
        
        estimate_no_residuals(svs=svs,rdir=rdir)
        print('Creating maps...')
        for sv in svs:
            _p(1,1,sv.name,flush=True)
            plot_map(sv,'-nr',rdir)
        write_csv(svs,os.path.join(rdir,'stats-nr.csv'))
        est_file = os.path.join(rdir,'estimates-nr.h5')
        res_file = os.path.join(rdir,'residuals-nr.h5')
        save_residuals(svs,est_file,res_file)

        
        svs,knns,loo_knns = estimate_loo_knn(svs=svs,k=20,estimates_file=est_file,rdir=rdir)
        print('Creating maps...')
        for sv in svs:
            _p(1,1,sv.name,flush=True)
            plot_map(sv,'-knn',rdir)
        write_csv(svs,os.path.join(rdir,'stats-knn.csv'))
        #preprocess(svs=svs,rdir=rdir+'-traces')
    except KeyboardInterrupt as ki:
        print('Interrupted by user')
    except Exception as e:
        print('Exception:',e)
        traceback.print_tb(e.__traceback__)
    #test_imu_delays()