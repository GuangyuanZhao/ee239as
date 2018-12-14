import numpy as np
from scipy import misc
from scipy.optimize import minimize
from PIL import Image
import pandas as pd
from estimate_position import poincare, az_zen_dist
# from libtiff import TIFF
import matplotlib.pyplot as plt
from PIL import Image
from simulate import oceanstokes, oceanaop
from polarization import *
from stats_utils import *
from scipy.optimize import minimize
# os.environ['PROJ_LIB'] = 'D:/NoSpace/Anaconda3/share/proj'
import os
import glob
import fnmatch

cwd = os.getcwd()


def _p(v, l, *args, **kw):
    """print if v, prefix l centered dots"""
    if v:
        if l: print('\xB7' * l, end='')  # B7 is centered dot
        print(*args, **kw)


def import_data(dir_path: object, ext: object) -> object:
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
    dic = file_dic
    i0 = np.asarray(dic["180"])
    i45 = np.asarray(dic["45"])
    i90 = np.asarray(dic["90"])
    i135 = np.asarray(dic["135"])
    return (i0, i45, i90, i135)


    ## Compute Stokes Parameters ##


def compute_stokes(i0, i45, i90, i135):
    #    pre-allocate an array to hold stokes vectors
    # s = np.ndarray(i0.shape + (3,))
    s = np.zeros((4,))
    # for inten in i0,i90,i45,i135
    #     for i in range (75,95)
    #     C{i - 74} = find(abs(inten - i) <= 2)
    #
    #     for i=1:length(C)
    #         bbb(i) = np.len(C{i})
    #
    #     index_shift = find(bbb == max(bbb))
    #     int_median = index_shift + 74
    #     average = mean(crop_new(C{index_shift}))
    if ndim([i0,i45,i90,i135])==1:
        s[..., 0] = np.average((i0 + i90 + i45 + i135) / 2.0)
        s[..., 1] = np.average(i0 - i90)
        s[..., 2] = np.average(i45 - i135)
        s[..., 3] = 0
    else:
        s[..., 0] = np.average((i0 + i90 + i45 + i135) / 2.0, axis=(0, 1))
        s[..., 1] = np.average(i0 - i90, axis=(0, 1))
        s[..., 2] = np.average(i45 - i135, axis=(0, 1))
        s[..., 3] = 0
    return s
    # For calculating errors in two forms.


def sun_pos_error_l2(sun_head, sun_zen, cam_aop, cam_head, cam_pitch, m2):
    sim_aop = oceanaop(sun_head, sun_zen, cam_head, cam_pitch, m2)
    d = angle_diff(cam_aop, sim_aop, pi)
    return sum(d ** 2)


def sun_pos_error_l1(sun_head, sun_zen, cam_aop, cam_head, cam_pitch, m2):
    sim_aop = oceanaop(sun_head, sun_zen, cam_head, cam_pitch, m2)
    d = angle_diff(cam_aop, sim_aop, pi)
    return sum(abs(d))


# Rename function
sun_pos_error = sun_pos_error_l1
residuals_name = 'residuals-L1.h5'


def fit_sun_head_to_aop(aop, head, pitch, ridx=1.33, verbose=False):
    # fit the sun heading to the data by assuming the sun zenith angle is pi/4:
    # TODO: the dimension of x0 should be the number of variables. So how many variables here?
    fit = minimize(sun_pos_error, x0=np.asarray([0]), args=(pi / 4, aop, head, pitch, ridx), bounds=[(-2 * pi, 2 * pi)],
                   options={'gtol': 1e-6})
    return fit.x


def fit_sun_to_aop(aop, head, pitch, sun_head_guess=None, ridx=1.33, verbose=False, vl=0):
    """fit_sun_to_aop, find sun position corresponding to aops at different headings & pitches, no time passing
    aop: array of AoPs, radians
    head: array of headings, radians, same shape as aop
    pitch: array of pitchs, radians, same shape as head
    ridx: refractive index of water
    """
    v = verbose
    _p(v, vl, 'Fitting sun ', end='')
    if sun_head_guess is None:
        _p(v, 0, 'heading to AoP... ', end='', flush=True)
        # fit just the heading first to get a good initial point:
        sun_head_guess = fit_sun_head_to_aop(aop, head, pitch, ridx)
    _p(v, 0, 'heading and zenith to AoP... ', end='', flush=True)
    # now do both heading & zenith
    minfun = lambda x, *args: sun_pos_error(x[0], x[1], *args)
    # Originally x0=(sun_head_guess,pi/4)
    # TODO: verify modification
    fit = minimize(minfun, x0=np.asarray([pi, pi/2]),method='L-BFGS-B', args=(aop, head, pitch, ridx),
                   bounds=[(-2 * pi, 2 * pi), (.01, pi / 2)], options={'gtol': 1e-6})
    sun_hz = fit.x
    _p(v, 0, 'DONE', flush=True)
    return sun_hz


if __name__ == '__main__':
    print('preprocessing...')
    b=[]
    # b = 'C:\\study\pol_GPS\\11_24\\13_10无太阳直射光\stored data'
    if b==[]:
        # f = open("C:\\study\\pol_GPS\\11_24\\13_10无太阳直射光\\stored data\\s.csv")
        # data = pd.read_csv(f)
        # print(data)
        # i0, i45, i90, i135 =data[0,0],data[1,0],data[2,0],data[3,0]
        # print (i0)
        # i0, i45, i90, i135 =[90.151519594495530,80.054397761443370,75.302827106820910,81.887167714389920] #indirect
        # i0, i45, i90, i135 =[1.008080338655501e+02,93.126594518204470,98.725996204933580,94.917654251316780]
        # i0, i45, i90, i135 =[93.126594518204470, 86.638399222294230, 94.917654251316780, 1.008080338655501e+02] #direct
        # i0, i45, i90, i135 =[1.008080338655501e+02,93.126594518204470,98.725996204933580,94.917654251316780] #direct2

        i0, i45, i90, i135 =[69.07496194999997,73.06496132000001,81.77103290000001,63.80920486] #weixi
    else:
        i0, i45, i90, i135 = import_data(b, 'tiff')
    s = compute_stokes(i0, i45, i90, i135)
    print('s=', s)
    # print(s/np.linalg.norm(s))
    (s0, p, a) = poincare(s)  # Caculating the intensity, degree of pol, angle of pol
    # sun_az, sun_zen, cam_head, cam_elev = math.radians(183), np.arccos(0.5924), math.radians(166), 0
    # sun_az, sun_zen, cam_head, cam_elev = math.radians(204.9), np.arccos(0.5182), math.radians(206), 0 # indirect illumination
    sun_az, sun_zen, cam_head, cam_elev = math.radians(201.04), np.arccos(0.5362), math.radians(200), math.radians(0) # indirect illumination
    stokes_theoretical = oceanstokes(sun_az, sun_zen, cam_head, cam_elev)
    print('s_theoretical=', np.squeeze(stokes_theoretical))
    aop_theoretical = Stokes.aop(np.squeeze(stokes_theoretical), axis=0)
    aop = Stokes.aop(s, axis=0)
    # aop=-0.5
    print('aop=', aop)
    print('aop_theoretical=', aop_theoretical)
    fit = fit_sun_to_aop(aop, cam_head, pitch=cam_elev, ridx=1.33, verbose=False)
    print(fit)
    print('Calculated',math.degrees(fit[0]),math.degrees(fit[1]))
    print('Theoretical',math.degrees(sun_az),math.degrees(sun_zen))
    # dis=az_zen_dist(np.deg2rad([214,89]),np.deg2rad([201,56]))
    # print(dis)
