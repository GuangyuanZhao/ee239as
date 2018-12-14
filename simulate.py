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
from polarization import *

def water_index(salinity,temperature,wavelen):
    """Index of refraction of water given salinity, temperature, and wavelength of light
    Based on "Empirical equation for the index of refraction of seawater" by Xiaohong Quan and Edward S. Fry, Applied Optics Vol. 34, No. 18, 1995
    parameters:
        salinity is in parts/1000, from 0 to 35
        temperature is in C, from 0 to 30
        wavelen is in nm, from 400 to 700
    returns:
        index of refraction
    """
    n0, n1, n2, n3, n4 = 1.31405, 1.779e-4, -1.05e-6, 1.6e-8, -2.02e-6
    n5, n6, n7, n8, n9 = 15.868,  0.01155,  -0.00423, -4382,  1.1455e6
    S = salinity
    T = temperature
    T2 = T**2
    L = wavelen
    return n0 + (n1 + n2*T + n3*T2)*S + n4*T2 + (n5 + n6*S+n7*T)/L + n8/(L**2) + n9/(L**3)

def p_wave_slope_ab(a,b,w,slick=False):
    """probability of a surface slope for wind speed w
    from "Measurement of the Roughness of the Sea Surface
    from Photographs of the Sun's Glitter" by Charles Cox and Walter Munk, 1954

    a is the azimuth angle of the steepest surface slope from the direction pointing into the wind
    b is the angle of the steepest slope
    w is the wind speed in m/s measured at 41 ft above sea level
    """
    zx = sin(a)*tan(b) #cross wind slope (right when looking up-wind)
    zy = cos(a)*tan(b) #up-wind slope
    return p_wave_slope(zx,zy,w,slick)

def p_wave_slope(zx,zy,w,slick=False):
    """probability of a surface slope for wind speed w
    from "Measurement of the Roughness of the Sea Surface
    from Photographs of the Sun's Glitter" by Charles Cox and Walter Munk, 1954

    zx is the slope in the cross-wind direction (right when looking up-wind)
    zy is the slope in the up-wind direction
    w is the wind speed in m/s measured at 41 ft above sea level
    """
    if slick:
        sc = sqrt(0.003 + 0.84e-3*w) # +/- 0.002
        su = sqrt(0.005 + 0.78e-3*w) # +/- 0.002
        c21 = 0.00 # +/- 0.02
        c03 = 0.02 # +/- 0.05
        c40 = 0.36 # +/- 0.24
        c22 = 0.10 # +/- 0.05
        c04 = 0.26 # +/- 0.31
    else:
        sc = sqrt(0.003 + 1.92e-3*w) # +/- 0.002
        su = sqrt(0.000 + 3.16e-3*w) # +/- 0.002
        c21 = 0.01 - 0.0086*w # +/- 0.03
        c03 = 0.04 - 0.033*w # +/- 0.12
        c40 = 0.40 # +/- 0.23
        c22 = 0.12 # +/- 0.06
        c04 = 0.23 # +/- 0.41
    xi = zx/sc
    et = zy/su
    xi2 = xi**2
    et2 = et**2
    return exp(-0.5*(xi2 + et2))*(1 - c21*(xi2-1)*et/2 - c03*(et**3 - 3*et)/6 + c40*(et2**2 - 6*et2 + 3)/24 + c22*(xi2-1)*(et2-1)/4 + c04*(et2**2-6*et2+3)/24 )/(2*pi*sc*su)

def Mxform(x1,y1,x2,y2):
    """Transform a stokes vector from coordinates x1,y1 to x2,y2"""
    return Jones.toMueller([[dot(x2,x1), dot(x2, y1)], [dot(y2,x1), dot(y2,y1)]])

def norm_cross(a,b):
    c = cross(a,b)
    return c/norm(c)

def Mrotv(k,x1,x2):
    """Rotation of a wave w/ propagation vector k rotating from x-axis x1 to x2"""
    y1 = norm_cross(k,x1)
    y2 = norm_cross(k,x2)
    return Mxform(x1,y1,x2,y2)

def vector_angle(a,b):
    return arccos(sum(a*b,-1)/sqrt(sum(a**2,-1)*sum(b**2,-1)))

def oceansim(sun_az,sun_zen,cam_head,cam_elev=0,m2=1.33,npart=1.08,mu=3.483, debug=True):
    """Compute Stokes vector observed by a camera facing "head" direction with sun azimuth and zenith
    :rtype: object
    """
    #npart is particle index for volume scattering function (real, > 1)
    #mu is Junge slope for volume scattering function (3-5)
    #in this system, the x is east, y is north, z is up
    #azimuth/heading angle is measured eastward from north
    #zenith angle is measured down from vertical
    #elevation angle is measured up from horizontal

    n = array([0,0,1]) # water surface normal vector
    m1 = 1.0 # air index of refraction
    #vector from sun:
    ki = -array([sin(sun_az)*sin(sun_zen), cos(sun_az)*sin(sun_zen), cos(sun_zen)])
    xi = norm_cross(n,ki)
    #transmitted sunlight
    #tx, ty are the transmission amplitude coefficients in the xt, yt directions
    kt,tx,ty = Fresnel.transmission(ki,n,m1,m2)
    xt = xi
    #vector to camera
    kc = -array([sin(cam_head)*cos(cam_elev),cos(cam_head)*cos(cam_elev),sin(cam_elev)])*norm(kt)
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
    b = Scattering.vspf_fournier(th_s,npart,mu)
    #  transform to camera's horizontal and up vectors
    mm4 = Mxform(xsc,ys, xc,yc)
    #Combined: mm4 . (b*mm3) . mm2 . mm1
    m = mm4.dot(b*mm3.dot(mm2.dot(mm1)))
    #stokes vector
    s = m.dot([1,0,0,0])
    #Caculating with considering the
    p=0
    if p==1:
        m_rot=Mueller.rotator(math.radians(10))
        s = m_rot.dot(m.dot([1,0,0,0]))
    if debug:
        return s,m,(ki,xi),(kt,xt,xst),(kc,xc,xsc),(mm1,mm2,mm3,b,mm4)
    else:
        return s,m

def oceanstokes(sun_az,sun_zen,cam_head,cam_elev=0,m2=1.33,npart=1.08,mu=3.483):
    """Compute stokes vectors"""
    b = broadcast(sun_az,sun_zen,cam_head,cam_elev,m2,npart,mu)
    st = empty(b.shape+(4,))
    st_flat = st.reshape((-1,4))
    for i,x in enumerate(b):
        st_flat[i] = oceansim(*x)[0]
    # here I change st to st_flat. I don't know the influence.
    return st_flat

def oceanaop(sun_az,sun_zen,cam_head,cam_elev=0,m2=1.33,npart=1.08,mu=3.483):
    """Compute aop"""
    return Stokes.aop(oceanstokes(sun_az,sun_zen,cam_head,cam_elev,m2,npart,mu),-1)
