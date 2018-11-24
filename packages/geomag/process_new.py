import numpy as np
from scipy import misc
from scipy.optimize import minimize
from PIL import Image
# from libtiff import TIFF
from matplotlib import *
from PIL import Image

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
#os.environ['PROJ_LIB'] = 'D:/NoSpace/Anaconda3/share/proj'
from mpl_toolkits.basemap import Basemap
from pint import UnitRegistry
import yaml
import os
import glob
import fnmatch

cwd = os.getcwd()

def _p(v,l,*args,**kw):
    """print if v, prefix l centered dots"""
    if v:
        if l: print('\xB7'*l,end='') #B7 is centered dot
        print(*args,**kw)


def import_data(dir_path, ext):
    root = dir_path
    file_dic = {}
    for file in glob.glob(os.path.join(root, "*." + ext)):
        file_path = os.path.join(root, file)
        if fnmatch.fnmatch(file, "*180*"):
            file_dic['180'] = Image.open(file_path)
            print("Found 180")
        elif fnmatch.fnmatch(file, "*45*"):
            file_dic['45'] = Image.open(file_path)
            print("Found 45")
        elif fnmatch.fnmatch(file, "*90*"):
            file_dic['90'] = Image.open(file_path)
            print("Found 90")
        elif fnmatch.fnmatch(file, "*135*"):
            file_dic['135'] = Image.open(file_path)
            print("Found 135")
        else:
            print("Do not us this file:" + file)

    return file_dic






    ## Compute Stokes Parameters ##
# split stack into components so that we can use them in numexpr
#     i0= Image.open('C:/Users/Guangyuan/Desktop/pol2','average_180.tiff')
#     i45= Image.open('C:/Users/Guangyuan/Desktop/pol2/average_45.tiff')
#     i90= Image.open('C:/Users/Guangyuan/Desktop/pol2/average_90.tiff')
#     i135= Image.open('C:/Users/Guangyuan/Desktop/pol2/average_135.tiff')

    #    pre-allocate an array to hold stokes vectors
    # s = np.ndarray(i0.shape + (3,))
    # s[..., 0] = (i0+i90+i45+i135)/2
    # s[..., 1] = i0-i90
    # s[..., 2] = i45-i135
    # print(s)




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


if __name__ == '__main__':
    print('Creating maps...')
    b='C:\\Users\\Guangyuan\\Desktop\\pol2'
    dic=import_data(b,'tiff')
    i0 = np.asarray(dic["180"])
    i45=np.asarray(dic["45"])
    i90=np.asarray(dic["90"])
    i90=np.asarray(dic["90"])


    print(i0)





