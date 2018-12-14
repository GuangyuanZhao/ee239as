import numpy as np
from scipy import optimize
from simulate import norm_cross, Mrotv, vector_angle
import math
from stats_utils import angle_diff
from packages.sunpos import sunpos
import datetime
import pytz
import os
import glob
from random import randint
import matplotlib.pyplot as plt
from skimage.color import rgb2gray
from PIL import Image
import pickle
import random
from geomag.emm import EMMMesh


# Hyperparameter
rfc_idx = 1.33
LATITUDE = 34.05
LONGITUTE = -118.24
# TIME = "2018-11-19 19:51:00"
# TIME = "2018-11-24 20:54:00"
# TIME = "2018-12-01 21:28:00"
# TIME = "2018-12-01 22:09:00"
TIME = "2018-12-13 19:21:00"
SIZE = 220
mu = 3.483 #3-5
npart = 1.08 #1-3
H_Min = 900
W_Min = 1583


def az_zen_dist(p0,p1):
    """angular distance between two points in azimuth,zenith format"""
    #formula comes from translating points into cartesian coordinates
    #taking the dot product to get the cosine between the two vectors
    #then arccos to return to angle, and simplify everything assuming real inputs
    a0,z0 = p0[0], p0[1]
    a1,z1 = p1[...,0], p1[...,1]
    return np.arccos(np.cos(z0)*np.cos(z1)+np.cos(a0-a1)*np.sin(z0)*np.sin(z1))

class Fresnel:
    """
    For all functions:
    mi,mt = indices of refraction of incident and transmit sides of interface
    ki,kt = wave-vectors, magnitude = 2*pi*m / wavelength
    n = normal vector of interface
    """
    @staticmethod
    def _vdot(a,b,axis=-1,keepdims=True):
        """sum-product along axis"""
        return np.sum(a*b,axis=axis,keepdims=keepdims)

    @staticmethod
    def _reflect(ki,n):
        #use inner instead of dot for broadcasting
        #k_a?
        #TODO: why times 2？
        return ki-2*n*Fresnel._vdot(n,ki)

    @staticmethod
    def reflect(ki,n):
        """Reflect a plane wave with propagation vector ki off an interface with unit normal vector n.
        The normal vector points towards the incident side of the surface.
        """
        #Transform into array
        ki,n = map(np.asarray, (ki,n))
        return Fresnel._reflect(ki,n)

    @staticmethod
    def _refract(ki,n,mi,mt):
        #use inner() instead of dot() for broadcasting
        #This is k_a
        kll = ki-n*Fresnel._vdot(ki,n)
        #This is k_t
        return kll - n*np.sqrt(Fresnel._vdot(ki,ki)*(mt/mi)**2 - Fresnel._vdot(kll,kll))
        #st = (mi/mt)*np.cross(np.cross(n,ki),n)
        #return st - n*np.sqrt(1-np.sum(st**2,axis=-1))

    @staticmethod
    def refract(ki,n,mi,mt):
        """Refract a plane wave with wave vector ki through an interface between mi and mt with unit normal vector n."""
        ki,n,mi,mt = map(np.asarray, (ki,n,mi,mt))
        return Fresnel._refract(ki,n,mi,mt)

    @staticmethod
    def _refl_trans(ki,kt,n,mi,mt):
        ss = Fresnel._vdot(ki + kt, n,keepdims=False)
        pp = Fresnel._vdot(mt**2*ki + mi**2*kt, n,keepdims=False)
        rs = Fresnel._vdot(ki - kt, n,keepdims=False) / ss
        rp = Fresnel._vdot(mt**2*ki - mi**2*kt, n,keepdims=False) / pp
        tt = 2*Fresnel._vdot(ki,n,keepdims=False)
        ts = tt / ss
        tp = mi*mt*tt / pp
        return rs, rp, ts, tp

    @staticmethod
    def _transmission(ki,kt,n,mi,mt):
        tt = 2*Fresnel._vdot(ki, n, keepdims=False)
        ts = tt / Fresnel._vdot(ki+kt, n, keepdims=False)
        tp = mi*mt*tt / Fresnel._vdot(mt**2*ki+mi**2*kt, n, keepdims=False)
        return ts, tp

    @staticmethod
    def transmission(ki,n,mi,mt):
        """Propagation vector and transmission coefficients (kt, ts, tp)"""
        ki,n,mi,mt = map(np.asarray,(ki,n,mi,mt))
        kt = Fresnel._refract(ki,n,mi,mt)
        return (kt,)+Fresnel._transmission(ki,kt,n,mi,mt)

    @staticmethod
    def _reflection(ki,kt,n,mi,mt):
        rs = Fresnel._vdot(ki - kt, n,keepdims=False) / Fresnel._vdot(ki + kt, n,keepdims=False)
        rp = Fresnel._vdot(mt**2*ki - mi**2*kt, n,keepdims=False) / Fresnel._vdot(mt**2*ki + mi**2*kt, n,keepdims=False)
        return rs, rp

    @staticmethod
    def reflection(ki,n,mi,mt):
        ki,n,mi,mt = map(np.asarray,(ki,n,mi,mt))
        kr = Fresnel._reflect(ki,n)
        kt = Fresnel._refract(ki,n,mi,mt)
        return (kr,)+Fresnel._reflection(ki,kt,n,mi,mt)

    @staticmethod
    def reflection_transmission(ki,n,mi,mt):
        ki,n,mi,mt = map(np.asarray,(ki,n,mi,mt))
        kr = Fresnel._reflect(ki,n)
        kt = Fresnel._refract(ki,n,mi,mt)
        rs,rp,ts,tp = Fresnel._refl_trans(ki,kt,n,mi,mt)
        return (kr,rs,rp,kt,ts,tp)

class Mueller:
    @staticmethod
    def rotator(angle):
        """
        Rotates Stokes Vector by angle radians from S1 to S2. This is the opposite of the Polarized Light book convention

        :return:
         _                                    _
        |   1  0             0             0   |
        |   0  cos(2*angle) -sin(2*angle)  0   |
        |   0  sin(2*angle)  cos(2*angle)  0   |
        |_  0  0             0             1  _|

        """

        c = np.cos(2*angle)
        s = np.sin(2*angle)
        return np.array([[1,0,0,0],[0,c,-s,0],[0,s,c,0],[0,0,0,1]])

    @staticmethod
    def rotate(mat,angle):
        """Rotate a Mueller matrix from S1 to S2 by angle radians"""
        return np.dot(Mueller.rotator(angle), np.dot(mat, Mueller.rotator(-angle)))

    @staticmethod
    def polarizer(px,py,angle=0):
        """
        :param px: t_s
        :param py: t_p
        :param angle:
        :return:
        corresponds to M_R in paper
         _                                       _
        |   px^2+py^2  px^2-py^2   0      0       |
        |   px^2-py^2  px^2+py^2   0      0       |
        |   0          0           2pxpy  0       |
        |_  0          0           0      2pxpy  _|

        """
        a = px**2+py**2
        b = px**2-py**2
        c = 2*px*py
        m = 0.5*np.array([[a,b,0,0],[b,a,0,0],[0,0,c,0],[0,0,0,c]])
        if angle != 0:
            return Mueller.rotate(m,angle)
        else:
            return m

    @staticmethod
    def retarder(phase,angle=0):
        """
        :param phase:
        :param angle:
        :return:
        m=
         _                                 _
        |   1  0   0           0            |
        |   0  1   0           0            |
        |   0  0   cos(angle)  sin(angle)   |
        |_  0  0  -sin(angle)  cos(angle)  _|

        """
        c = np.cos(phase)
        s = np.sin(phase)
        m = np.array([[1,0,0,0],[0,1,0,0],[0,0,c,s],[0,0,-s,c]])
        if angle != 0:
            return Mueller.rotate(m,angle)
        else:
            return m

    @staticmethod
    def depolarizer(dp):
        """
        A scaling of dimension 2,3,4
        """
        return np.array([[1,0,0,0],[0,dp,0,0],[0,0,dp,0],[0,0,0,dp]])

    @staticmethod
    def rayleigh(th,r,wl,a,n1,n2):
        """Rayleigh scattering Mueller matrix

        Paramters:
          th is scattering angle, radians
          r is the observation distance, m
          wl is wavelength in medium, m
          a is sphere radius, m
          n1 is medium (real) index of refraction
          n2 is sphere (complex) index of refraction
        Returns:
          The 4x4 real Mueller matrix for Rayleigh scattering
        Notes:
          The incident and scattered light waves are not in the same coordinate system
          If the incident wave is propagating along z in x,y,z = (1,0,0),(0,1,0),(0,0,1)
          Then the scattered wave will be propagating along z' in x',y',z' where
          x' = (cos(th),0,-sin(th)), y' = (0,1,0), z' = (sin(th),0,cos(th)).
          N.B. that S1 = 1 means polarized along +/- x, S1 = -1 means polarized along +/- y

          The equation for this matrix is given by:
            (k**4) * (a**6) * abs(n**2 - 1)**2/abs(n**2 + 2)**2 /(2 * r**2) * [[1 + c**2, -s**2, 0, 0],[-s**2, 1 + c**2, 0, 0],[0, 0, 2*c, 0],[0, 0, 0, 2*c]]
            Where k = 2*pi/wl, and n**2 = n2**2/n1**2, c = cos(th), s = sin(th)
        """
        c = np.cos(th)
        c2,s2 = c**2, np.sin(th)**2
        k = 2*np.pi/wl
        n_2 = n2**2/n1**2
        m = (k**4)*(a**6)*(abs(n_2-1)**2) / ((abs(n_2+2)**2) * 2 * (r**2))
        return m*np.array([[1+c2 , -s2  , 0   , 0],
                           [-s2  , 1+c2 , 0   , 0],
                           [0    , 0    , 2*c , 0],
                           [0    , 0    , 0   , 2*c]])

    @staticmethod
    def rayleigh_noscale(th):
        """Rayleigh scattering Mueller matrix, no physical scaling
        Paramters:
          th is scattering angle, radians
        Returns:
          The 4x4 real Mueller matrix for Rayleigh scattering, with no physical scaling
        Notes:
          The incident and scattered light waves are not in the same coordinate system
          If the incident wave is propagating along z in x,y,z = (1,0,0),(0,1,0),(0,0,1)
          Then the scattered wave will be propagating along z' in x',y',z' where
          x' = (cos(th),0,-sin(th)), y' = (0,1,0), z' = (sin(th),0,cos(th)).
          N.B. that S1 = 1 means polarized along +/- x, S1 = -1 means polarized along +/- y

          The equation for this matrix is given by:
            0.5*[[c**2+1, c**2-1, 0, 0],[c**2-1, c**2+1, 0, 0],[0, 0, 2*c, 0],[0, 0, 0, 2*c]]
            Where c = cos(th)
        """
        c = np.cos(th)
        c2 = c**2
        return 0.5*np.array([[c2+1 , c2-1  , 0   , 0],
                           [c2-1  , c2+1 , 0   , 0],
                           [0    , 0    , 2*c , 0],
                           [0    , 0    , 0   , 2*c]])

    @staticmethod
    def rayleigh_norm(th):
        """Rayleigh scattering Mueller matrix, normalized by intensity
        Paramters:
          th is scattering angle, radians
        Returns:
          The 4x4 real, normalized Mueller matrix for Rayleigh scattering
        Notes:
          The incident and scattered light waves are not in the same coordinate system
          If the incident wave is propagating along z in x,y,z = (1,0,0),(0,1,0),(0,0,1)
          Then the scattered wave will be propagating along z' in x',y',z' where
          x' = (cos(th),0,-sin(th)), y' = (0,1,0), z' = (sin(th),0,cos(th)).
          N.B. that S1 = 1 means polarized along +/- x, S1 = -1 means polarized along +/- y

          The equation for this matrix is given by:
            [[1, (c**2-1)/(c**2+1), 0, 0],[(c**2-1)/(c**2+1), 1, 0, 0],[0, 0, 2*c/(c**2+1), 0],[0, 0, 0, 2*c/(c**2+1)]]
            Where c = cos(th), s = sin(th)
        """
        c = np.cos(th)
        c2 = c**2
        a = (c2-1)/(c2+1)
        b = 2*c/(c2+1)
        return np.array([[1,a,0,0],[a,1,0,0],[0,0,b,0],[0,0,0,b]])



    _pauli = np.array([[[1,0],[0,1]],[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]])
    #@staticmethod
    #def valid(mat):
    #    pass

class Scattering:
    @staticmethod
    def vspf_fournier(th,n,mu):
        """Volume scattering phase function, Fournier & Forand (1999), Fournier & Jonasz (1999)
        Parameters:
            th = scattering angle
            n = index of refraction of particle, only the real part is used
            mu = slope of hyperbolic distribution of particle sizes, typically 3-5
        """
        n = np.real(n)
        d = 4*np.sin(th/2)**2 / (3*(n-1)**2)
        d_180 = 4*np.sin(np.pi/2) / (3*(n-1)**2)
        v = (3-mu)/2
        dv = d**v
        d_180v = d_180**v
        d1 = 1-d
        dv1 = 1-dv

        a = 1/(4*np.pi*dv*d1**2)
        b = v*d1-dv1+(d*dv1-v*d1)*np.sin(th/2)**(-2)
        c = (1-d_180v)*(3*np.cos(th)**2 - 1)/(16*np.pi*(d_180-1)*d_180v)
        return a*b+c

    @staticmethod
    def bsf_fournier(n,mu):
        """Backscattering fraction
        Parameters:
            n = index of refraction of particles, only the real part is used
            mu = slope of hyperbolic distribution of particle sizes, typically 3-5
        """
        n = np.real(n)
        d_90 = 4*np.sin(np.pi/4)/(3*(n-1)**2)
        v = (3-mu)/2
        d_90v = d_90**v
        a = 1 - d_90**(v+1) - 0.5*(1-d_90v)
        b = (1-d_90)*d_90v
        return 1 - a/b

class Jones:
    @staticmethod
    def toMueller(mat):
        mat = np.array(mat,copy=False)
        A = np.array([[1,0,0,1],[1,0,0,-1],[0,1,1,0],[0,1j,-1j,0]])
        return 0.5*np.real(A.dot(np.kron(mat, mat.conj())).dot(A.T.conj()))

    @staticmethod
    def toStokes(vec):
        vec = np.array(vec,copy=False)
        A = np.array([[1,0,0,1],[1,0,0,-1],[0,1,1,0],[0,1j,-1j,0]])
        return np.real(A.dot(np.kron(vec,vec.conj())))

    @staticmethod
    def rotate(mat,angle):
        """Rotate a Jones matrix from x to y by angle radians"""
        return np.dot(Jones.rotator(angle), np.dot(mat, Jones.rotator(-angle)))

    @staticmethod
    def rotator(angle):
        """Rotates electric field by angle radians from x to y. This is the opposite of the Polarized Light book convention"""
        c = np.cos(angle)
        s = np.sin(angle)
        return np.array([[c,-s],[s,c]])

    @staticmethod
    def polarizer(px,py,angle=0):
        """Polarizer/diattenuator with px, py as x, y amplitude coefficients and rotation angle (radians)"""
        M = np.array([[px,0],[0,py]])
        if angle != 0:
            return Jones.rotate(M,angle)
        else:
            return M

    @staticmethod
    def retarder(phase,angle=0):
        """Retarder with phase and rotation angle both in radians"""
        r = np.exp(1j*phase/2)
        R = np.array([[r,0],[0,np.conj(r)]])
        if angle != 0:
            return Jones.rotate(R,angle)
        else:
            return R

    @staticmethod
    def rayleigh(th,r,wl,a,n1,n2):
        """Rayleigh scattering jones matrix
        Paramters:
          th is scattering angle, radians
          r is the observation distance, m
          wl is wavelength in medium, m
          a is sphere radius, m
          n1 is medium (real) index of refraction
          n2 is sphere (complex) index of refraction
        Returns:
          The 2x2 complex Jones matrix for Rayleigh scattering
        Notes:
          The incident and scattered light waves are not in the same coordinate system
          If the incident wave is propagating along z in x,y,z = (1,0,0),(0,1,0),(0,0,1)
          Then the scattered wave will be propagating along z' in x',y',z' where
          x' = (cos(th),0,-sin(th)), y' = (0,1,0), z' = (sin(th),0,cos(th)).

          The equation for this matrix is given by:
            (k**2) * (a**3) * (n**2 - 1)/(n**2 + 2) * 1/r * [[cos(th), 0],[0, 1]]
            Where k = 2*pi/wl, and n**2 = n2**2/n1**2
        """
        k = 2*np.pi/wl
        n_2 = n2**2/n1**2
        return ((k**2)*(a**3)*((n_2-1)/(n_2+2))/r)*np.array([[np.cos(th), 0],[0,1]])

class Stokes:
    @staticmethod
    def intensity(s,axis=0):
        return np.take(s,0,axis)
    @staticmethod
    def aop(s,axis=0):
        """
        Angle of polarizatioin
        :param s:
        :param axis:
        :return:
        """
        s1 = np.take(s,1,axis)
        s2 = np.take(s,2,axis)
        return 0.5*np.arctan2(s2,s1)
    @staticmethod
    def dop(s,axis=0):
        """
        Degree of polarization
        :param s:
        :param axis:
        :return:
        """
        s0 = np.take(s,0,axis)
        s123 = np.take(s,(1,2,3),axis)
        return np.linalg.norm(s123,axis=axis)/s0
    @staticmethod
    def dolp(s,axis=0):
        """
        Degree of linear polarization
        :param s:
        :param axis:
        :return:
        """
        s0 = np.take(s,0,axis)
        s12 = np.take(s,(1,2),axis)
        return np.linalg.norm(s12,axis=axis)/s0
    @staticmethod
    def docp(s,axis=0):
        """
        degree of circular polarization
        :param s:
        :param axis:
        :return:
        """
        s0 = np.take(s,0,axis)
        s3 = np.take(s,3,axis)
        return abs(s3)/s0
    @staticmethod
    def ella(s,axis=0):
        s0 = np.take(s,0,axis)
        s3 = np.take(s,3,axis)
        return 0.5*np.arcsin(s3/s0)
    @staticmethod
    def to_poincare(s,axis=0):
        """intensity, DoP, AoP, Ellipticity"""
        d = list(np.shape(s))
        d[axis] = 4
        p = np.empty(shape=d)
        pv = p.swapaxes(0,axis)
        pv[0] = Stokes.intensity(s,axis)
        pv[1] = Stokes.dop(s,axis)
        pv[2] = Stokes.aop(s,axis)
        pv[3] = Stokes.ella(s,axis)
        return p
    @staticmethod
    def to_poincare_linear(s,axis=0):
        """intensity, DoLP, AoP"""
        d = list(np.shape(s))
        d[axis] = 3
        p = np.empty(shape=d)
        pv=p.swapaxes(0,axis)
        pv[0] = Stokes.intensity(s,axis)
        pv[1] = Stokes.dolp(s,axis)
        pv[2] = Stokes.aop(s,axis)
        return p


def Mxform(x1,y1,x2,y2):
    """Transform a stokes vector from coordinates x1,y1 to x2,y2"""
    return Jones.toMueller([[np.dot(x2,x1), np.dot(x2, y1)], [np.dot(y2,x1), np.dot(y2,y1)]])

#Functions that compute Stokes vector given position of sun
def oceansim(sun_az,sun_zen,cam_head,cam_elev=0,m2=1.33,npart=1.08,mu=3.483, debug=True):
    """
    Compute Stokes vector observed by a camera facing "head" direction with sun azimuth and zenith
    -in this system, the x is east, y is north, z is up
    -azimuth/heading angle is measured eastward from north
    -zenith angle is measured down from vertical
    -elevation angle is measured up from horizontal
    Args:
        sun_az: azimuth
        sun_zen: zenith
        cam_head: heading angle?
        cam_elev: elevation
        m2: index of refraction of water
        npart: particle index for volume scattering function (real, > 1)
        mu: Junge slope for volume scattering function (3-5)
        debug:
    Return:
        s: stokes vector
    """

    #Water surface norm
    n = np.array([0,0,1])
    m1 = 1.0
    #vector from sun:
    ki = -np.asarray([np.sin(sun_az)*np.sin(sun_zen),
                      np.cos(sun_az)*np.sin(sun_zen),
                      np.cos(sun_zen)])
    xi = norm_cross(n,ki)
    #transmitted sunlight
    #tx, ty are the transmission amplitude coefficients in the xt, yt directions
    kt,tx,ty = Fresnel.transmission(ki,n,m1,m2)
    xt = xi
    #vector to camera
    kc = -np.asarray([np.sin(cam_head)*np.cos(cam_elev),
                      np.cos(cam_head)*np.cos(cam_elev),
                      np.sin(cam_elev)])*np.linalg.norm(kt)
    xc = norm_cross(n, kc) #right
    yc = norm_cross(kc, xc) #up
    #vectors for scattering
    ys = norm_cross(kt, kc) # y-axis of scattering event
    xst = norm_cross(ys, kt) # x-axis of scattering event relative to transmitted sunlight
    xsc = norm_cross(ys, kc) # x-axis of scattering event relative to camera
    #Mueller matrices
    #  transmission through water surface:
    mm1 = Mueller.polarizer(tx,ty)
    #  rotate to scattering plane
    mm2 = Mrotv(kt,xt,xst)
    #  scatter
    th_s = vector_angle(kt,kc)
    #mm3 = Mocean(rad2deg(th_s)) #using Empirical ocean scattering
    mm3 = Mueller.rayleigh_norm(th_s) #normalized Rayleigh scattering matrix
    #b = Scattering.bsf_fournier(npart,mu)
    b = Scattering.vspf_fournier(th_s,npart,mu)
    #  transform to camera's horizontal and up vectors
    mm4 = Mxform(xsc,ys, xc,yc)
    #Combined: mm4 . (b*mm3) . mm2 . mm1
    m = mm4.dot(b*mm3.dot(mm2.dot(mm1)))
    #stokes vector
    s = m.dot([1,0,0,0])
    if debug:
        return s,m,(ki,xi),(kt,xt,xst),(kc,xc,xsc),(mm1,mm2,mm3,b,mm4)
    else:
        return s,m

def oceanstokes(sun_az,sun_zen,cam_head,cam_elev=0,m2=1.33,npart=1.08,mu=3.483):
    """
    Compute stokes vectors
    Args:
        sun_az: azimuth scalar
        sun_zen: zenith scalar
        cam_head: heading in radian
        m2: refractive index of water
        npart: particle index
    Returns:
        st: stokes vector
    """

    #TODO: ???
    b = np.broadcast(sun_az,sun_zen,cam_head,cam_elev,m2,npart,mu)
    st = np.empty(b.shape+(4,))
    st_flat = st.reshape((-1,4))
    for i,x in enumerate(b):
        s, m, incident, refraction, camera, para = oceansim(*x, debug=True)
        st_flat[i] = s
    return st

def oceanaop(sun_az,sun_zen,cam_head,cam_elev=0,m2=1.33,npart=1.08,mu=3.483):
    """Compute aop"""
    stokes = oceanstokes(sun_az,sun_zen,cam_head,cam_elev,m2,npart,mu)
    aop = Stokes.aop(stokes, -1)
    return aop

#Estimate sun position based on model
def sun_pos_error_l2(sun_head,sun_zen,cam_aop,cam_head,cam_pitch,m2):
    sim_aop = oceanaop(sun_head,sun_zen,cam_head,cam_pitch,m2)
    d = angle_diff(cam_aop,sim_aop,np.pi)
    return np.sum(d**2)

def sun_pos_error_l1(sun_head,sun_zen,cam_aop,cam_head,cam_pitch,m2):
    sim_aop = oceanaop(sun_head,sun_zen,cam_head,cam_pitch,m2)
    d = angle_diff(cam_aop,sim_aop,np.pi)
    return np.sum(abs(d))

sun_pos_error = sun_pos_error_l2

def fit_sun_head_to_aop(aop,head,pitch,ridx=1.33,verbose=False):
    #fit the sun heading to the data by assuming the sun zenith angle is pi/4:
    #TODO: the dimension of x0 should be the number of variables. So how many variables here?
    fit = optimize.minimize(sun_pos_error,x0=np.asarray([np.pi]),args=(np.pi/4,aop,head,pitch,ridx),bounds=[(-2*np.pi,2*np.pi)],options={'gtol':1e-6})
    return fit.x

def fit_sun_to_aop(aop,head,pitch,sun_head_guess=None,ridx=1.33,verbose=False,vl=0):
    """fit_sun_to_aop, find sun position corresponding to aops at different headings & pitches, no time passing
    aop: array of AoPs, radians
    head: array of headings, radians, same shape as aop
    pitch: array of pitchs, radians, same shape as head
    ridx: refractive index of water
    """
    v = verbose
    if sun_head_guess is None:
        #fit just the heading first to get a good initial point:
        sun_head_guess = fit_sun_head_to_aop(aop,head,pitch,ridx)
        pass
    #now do both heading & zenith
    minfun = lambda x,*args: sun_pos_error(x[0],x[1],*args)
    #Originally x0=(sun_head_guess,pi/4)
    #TODO: verify modification
    fit = optimize.minimize(minfun,x0=np.asarray([np.pi,np.pi/4]),args=(aop,head,pitch,ridx),bounds=[(0,2*np.pi),(.01,np.pi)],options={'gtol':1e-6})
    sun_hz = fit.x
    # print('DONE')
    return sun_hz

def greedy_minimization(aop, head, pitch, sun_head_guess=None, ridx=1.33):
    print("Start minimization!")
    min_err = None
    for az in range(1,360):
        for zen in range(1,90):
            err = sun_pos_error(np.deg2rad(az), np.deg2rad(zen), aop, head, pitch, ridx)
            if min_err is None:
                min_err = err
                opt_value = (az, zen)
            elif min_err<err:
                min_err = err
                opt_value = (az, zen)
            elif min_err==err:
                print("SAME ERROR AT (%f,%f) AND (%f,%f)"%(opt_value[0], opt_value[1], az, zen))
                print("ERROR: %f" % min_err)
            else:
                pass
    return opt_value

def declination(lat,lon,elev,t,gm,radians=True):
    b = np.broadcast(lat, lon, elev, t)
    res = np.empty(b.shape)
    for i,x in enumerate(b):
        la,lo,el,tt = x
        res.flat[i] = gm.declination(la,lo,el,gm.decimal_year(tt))
    if radians:
        res = np.deg2rad(res)
    return res

def sunpos_mag(t,lat,lon,elev,gm,temp=None,press=None,radians=True):
    """Observed sun heading, zenith -- using magnetic N"""
    #az_zen is a (...,5) dimension ndarray
    az_zen = sunpos(t,lat,lon,elev,temp,press,radians=radians)
    decl = declination(lat,lon,elev,t,gm,radians)
    az_zen[...,0] -= decl
    #subtract declination to go from true N to magnetic N
    return az_zen

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
    aop = oceanaop(sun_head,sun_zen,head,pitch,1.33)
    return sun_zen,sun_head,aop

# def fit_sun_alt(aop,head,pitch,sun_head_guess=None,ridx=1.33,verbose=False,vl=0):
#     sun_az = Variable(1)
#     sun_zen = Variable(1)
#
#     objective = Minimize(sun_pos_error)
#     constraints = [0 <= sun_az <= 2 * np.pi,
#                    0 <= sun_zen <= np.pi / 2]
#     prob = Problem(objective, constraints)
#     result = prob.solve()
#     print(sun_az.value, sun_zen.value)
#     print('DONE\n')
#     return (sun_az.value, sun_zen.value)


#Variable: 2 scalar, sun_head, sun_zen
#function:

class SpinImage:
    def __init__(self, prefix, cam_head, cam_elev=0, temp=22, pres=1013.0, dt=0.0):
        self.lat = LATITUDE
        self.lon = LONGITUTE
        self.time = datetime.datetime.strptime(TIME, '%Y-%m-%d %H:%M:%S')
        self.elev = cam_elev
        self.head = cam_head
        self.temp = temp
        self.pressure = pres
        self.dt = dt
        self.aop = None
        self.dop = None
        self.prefix = prefix
        self.stokes = None
        self.folders = None
        self.model_stokes = None

    def get_dir(self):
        self.folders = {}
        for root, dirs, files in os.walk(self.prefix):
            for name in dirs:
                angle = int(name.split('_')[-1])
                self.folders[angle] = self.prefix+'/'+name
        return self.folders

    def read_images(self, denoise=True, read_all=False, i=1):
        if self.folders is None:
            self.get_dir()
        img_list = []
        for angle in sorted(self.folders.keys()):
            files = glob.glob(self.folders[angle] + '/*.jpg')
            img = []

            if not read_all:
                # Randomly read one image
                i_list = random.sample(range(len(files)), i)
                for j in i_list:
                    index = j
                    single_img = plt.imread(files[index])
                    single_img = np.asarray(0.299 * single_img[..., 0] + 0.587 * single_img[..., 1] + 0.114 * single_img[..., 2])
                    crop = single_img[H_Min:H_Min+SIZE, W_Min:W_Min+SIZE]
                    img.append(crop)
            else:
                # Read all 30 images
                for index in range(len(files)):
                    single_img = plt.imread(files[index])
                    single_img = np.asarray(0.299 * single_img[..., 0] + 0.587 * single_img[..., 1] + 0.114 * single_img[..., 2])
                    crop = single_img[H_Min:H_Min + SIZE, W_Min:W_Min + SIZE]
                    img.append(crop)

            # Calculate average over 30 images and 220x220 pixels. img is a scalar.
            img = np.average(img, axis=0)
            img_list.append(img)

        # Sort image list
        # img_list = [img_list[i] for i in [2,3,0,1]]

        if denoise:
            img_list = self._denoise(img_list)
        return img_list

    def _denoise(self, img_list):
        img_list_denoised = []
        for i in range(4):
            img = img_list[i]
            mask = []
            for j in range(int(np.min(img)), int(np.max(img))):
                mask.append(np.abs(img-j)<=3)

            k_argmax = np.argmax([np.sum(mask[k]) for k in range(len(mask))])
            average = np.sum(mask[k_argmax]*img)/np.sum(mask[k_argmax])
            img_list_denoised.append(average)
        return img_list_denoised

    def compute_stokes(self, im3g_list):
        self.stokes = np.ndarray(img_list[0].shape+(3,))
        self.stokes[...,0] = sum(img_list[i] for i in range(4))/4.0
        self.stokes[...,1] = img_list[0] - img_list[2]
        self.stokes[...,2] = img_list[1] - img_list[3]

    def compute_aop(self):
        if self.stokes is None:
            print('Calculating stokes vector first')
            img_list = self.read_images()
            self.compute_stokes(img_list)
        s1 = np.take(self.stokes, 1, axis=-1)
        s2 = np.take(self.stokes, 2, axis=-1)
        self.aop = 0.5*np.arctan2(s2,s1)

    def compute_dop(self):
        if self.stokes is None:
            print('Calculating stokes vector first...')
            img_list = self.read_images()
            self.compute_stokes(img_list)
        s0 = np.take(self.stokes, 0, axis=-1)
        s12 = np.take(self.stokes, (1,2), axis=-1)
        self.dolp = np.linalg.norm(s12, axis=-1)/s0
        pass

    def get_ground_truth(self, rad=True):
        az, zen, _, _, _ = sunpos(self.time, self.lat, self.lon,
                            self.elev, self.temp, self.pressure,
                            self.dt, rad)
        self.sun = np.asarray([az, zen])
        self.model_stokes = oceanstokes(self.sun[0], self.sun[1], self.head,
                    0, rfc_idx, npart=npart, mu=mu)
        s1 = np.take(self.model_stokes, 1, axis=-1)
        s2 = np.take(self.model_stokes, 2, axis=-1)
        self.model_aop = 0.5*np.arctan2(s2, s1)
        del s1, s2
        pass

    def compute_residual(self):
        if self.aop is None:
            self.compute_aop()
        if self.model_stokes is None or self.sun is None:
            self.get_ground_truth()
        self.sun_head_est = fit_sun_head_to_aop(self.aop, self.head, 0)
        self.cam_head_rel = angle_diff(self.head, self.sun_head_est)
        self.model_aop_res = angle_diff(self.aop, self.model_aop, np.pi)

    def estimate_sun_pos(self):
        self.compute_residual()
        # fit_sun_pos = fit_sun_to_aop(self.aop, self.head, 0,
        #                                   self.sun_head_est, rfc_idx)
        fit_sun_pos = fit_sun_to_aop(self.aop, self.head, 0,
                                          self.sun_head_est, rfc_idx)
        fit_sun_err = az_zen_dist(self.sun, np.asarray(fit_sun_pos))
        self.fit_sun_pos = fit_sun_pos
        self.fit_sun_err = fit_sun_err
        return fit_sun_pos, fit_sun_err


if __name__ == '__main__':
    # prefix = './data/Pol_11.47'
    # prefix = './data/11_24/12_54'
    # prefix = './data/12_01/group4'
    # prefix = './data/12_01/group6'
    prefix = "./data/12_13"

    print('Adding images together...')
    data1 = SpinImage(prefix, math.radians(167), 85)
    img_list = data1.read_images(denoise=True, read_all=True)
    data1.compute_stokes(img_list)
    sun_pos, error = data1.estimate_sun_pos()
    print("Estimated sun azimuth: %f, zenith: %f\n" % (np.rad2deg(sun_pos[0]), np.rad2deg(sun_pos[1])))
    print("Sun position error: %f" % np.rad2deg(error))

    # batch_num = 30
    # avg_dis = [0, 0]
    # avg_err = 0
    # for i in range(batch_num):
    #     data = SpinImage(prefix, math.radians(200), 85)
    #     img_list = data.read_images(denoise=True, read_all=False, i=5)
    #     data.compute_stokes(img_list)
    #     sun_pos, error = data.estimate_sun_pos()
    #     if i == 0:
    #         print("Ground truth: (%f deg, %f deg)" % (np.rad2deg(data.sun[0]), np.rad2deg(data.sun[1])))
    #     print("Batch %d, Estimated sun azimuth and zenith: (%f deg, %f deg),  Error: %f deg" % (i+1, np.rad2deg(sun_pos[0]), np.rad2deg(sun_pos[1]), np.rad2deg(error)))
    #     avg_dis[0] += np.rad2deg(abs(data.sun[0]-sun_pos[0]))
    #     avg_dis[1] += np.rad2deg(abs(data.sun[1]-sun_pos[1]))
    #     avg_err += np.rad2deg(error)
    #     del data
    # avg_dis[0] /= batch_num
    # avg_dis[1] /= batch_num
    # avg_err /= batch_num
    # print("Average Azimuth error: %f deg, Average Zenith error: %f deg, Average error: %f deg" % (avg_dis[0], avg_dis[1], avg_err))

    # print(np.rad2deg(sun_pos), np.rad2deg(error))

    # sun_az, sun_zen, cam_head, cam_elev = (math.radians(183), np.arccos(0.5924), math.radians(166), 0)
    # stokes = oceanstokes(sun_az, sun_zen, cam_head, cam_elev)
    # aop = oceanaop(sun_az, sun_zen, cam_head, cam_elev)