# The MIT License (MIT)
# 
# Copyright (c) 2016 Samuel Bear Powell
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np
cimport numpy as np
from libc.math cimport cos,sin,tan,sqrt,atan,atan2,asin,acos
from cpython cimport array
from datetime import datetime

cdef double deg2rad(double x):
    return x*np.pi/180.0

cdef double rad2deg(double x):
    return x*180.0/np.pi

cdef double polyval(p,double x):
    cdef double y = 0
    for i,v in enumerate(p):
        y *= x
        y += v
    return y

cdef tuple calendar_time(dt):
    try:
        x = dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second, dt.microsecond
        return x
    except AttributeError:
        try:
            return calendar_time(datetime.utcfromtimestamp(dt)) #will raise OSError if dt is not acceptable
        except:
            raise TypeError('dt must be datetime object or POSIX timestamp')

cdef double cjulian_day(dt):
    """Calculate the Julian Day from a datetime.datetime object in UTC"""
    # year and month numbers
    cdef int yr, mo, idy, hr, mn, sc, us
    yr, mo, idy, hr, mn, sc, us = calendar_time(dt)
    if mo <= 2:  # From paper: "if M = 1 or 2, then Y = Y - 1 and M = M + 12"
        mo += 12
        yr -= 1
    # day of the month with decimal time
    cdef double dy = idy + hr/24.0 + mn/(24.0*60.0) + sc/(24.0*60.0*60.0) + us/(24.0*60.0*60.0*1e6)
    # b is equal to 0 for the julian calendar and is equal to (2- A +
    # INT(A/4)), A = INT(Y/100), for the gregorian calendar
    cdef int a = int(yr / 100.0)
    cdef int b = 2 - a + int(a / 4.0)
    cdef double jd = int(365.25 * (yr + 4716)) + int(30.6001 * (mo + 1)) + dy + b - 1524.5
    return jd

cdef double julian_ephemeris_day(double jd,double deltat):
    """Calculate the Julian Ephemeris Day from the Julian Day and delta-time = (terrestrial time - universal time) in seconds"""
    return jd + deltat / 86400.0

cdef double julian_century(double jd):
    """Caluclate the Julian Century from Julian Day or Julian Ephemeris Day"""
    return (jd - 2451545.0) / 36525.0

cdef double julian_millennium(double jc):
    """Calculate the Julian Millennium from Julian Ephemeris Century"""
    return jc / 10.0

# Earth Periodic Terms
# Earth Heliocentric Longitude coefficients (L0, L1, L2, L3, L4, and L5 in paper)
cdef int _L5L = 1
cdef double* _L5A = [1]
cdef double* _L5B = [3.14]
cdef double* _L5C = [0.0]

cdef int _L4L = 3
cdef double* _L4A = [114,8,1]
cdef double* _L4B = [3.142,4.13,3.84]
cdef double* _L4C = [0.0,6283.08,12566.15]

cdef int _L3L = 7
cdef double* _L3A = [289, 35, 17, 3, 1, 1, 1]
cdef double* _L3B = [5.844, 0.0, 5.49, 5.2, 4.72, 5.3, 5.97]
cdef double* _L3C = [6283.076, 0.0, 12566.15, 155.42, 3.52, 18849.23, 242.73]

cdef int _L2L = 20
cdef double* _L2A = [52919, 8720, 309, 27, 16, 16, 10, 9, 7, 5, 4, 4, 3, 3, 3, 3, 3, 3, 2, 2]
cdef double* _L2B = [0.0,  1.0721,  0.867,  0.05,  5.19,  3.68,  0.76,  2.06,  0.83,  4.66,  1.03,  3.44,  5.14,  6.05,  1.19,  6.12,  0.31,  2.28,  4.38,  3.75]
cdef double* _L2C = [0.0,  6283.0758,  12566.152,  3.52,  26.3,  155.42,  18849.23,  77713.77,  775.52,  1577.34,  7.11,  5573.14,  796.3,  5507.55,  242.73,  529.69,  398.15,  553.57,  5223.69,  0.98]

cdef int _L1L = 34
cdef double* _L1A = [628331966747,  206059,  4303,  425,  119,  109,  93,  72,  68,  67,  59,  56,  45,  36,  29,  21,  19,  19,  17,  16,  16,  15,  12,  12,  12,  12,  11,  10,  10,  9,  9,  8,  6,  6]
cdef double* _L1B = [0.0,  2.678235,  2.6351,  1.59,  5.796,  2.966,  2.59,  1.14,  1.87,  4.41,  2.89,  2.17,  0.4,  0.47,  2.65,  5.34,  1.85,  4.97,  2.99,  0.03,  1.43,  1.21,  2.83,  3.26,  5.27,  2.08,  0.77,  1.3,  4.24,  2.7,  5.64,  5.3,  2.65,  4.67]
cdef double* _L1C = [0.0,  6283.07585,  12566.1517,  3.523,  26.298,  1577.344,  18849.23,  529.69,  398.15,  5507.55,  5223.69,  155.42,  796.3,  775.52,  7.11,  0.98,  5486.78,  213.3,  6275.96,  2544.31,  2146.17,  10977.08,  1748.02,  5088.63,  1194.45,  4694,  553.57,  3286.6,  1349.87,  242.73,  951.72,  2352.87,  9437.76,  4690.48]

cdef int _L0L = 64
cdef double* _L0A = [175347046,  3341656,  34894,  3497,  3418,  3136,  2676,  2343,  1324,    1273,  1199,  990,  902,  857,  780,  753,  505,  492,  357,  317,  284,    271,  243,  206,  205,  202,  156,  132,  126,  115,  103,  102,  102,  99,    98,  86,  85,  85,  80,  79,  71,  74,  74,  70,  62,  61,  57,  56,  56,  52,    52,  51,  49,  41,  41,  39,  37,  37,  36,  36,  33,  30,  30,  25]
cdef double* _L0B = [0.0,  4.6692568,  4.6261,  2.7441,  2.8289,  3.6277,  4.4181,  6.1352,  0.7425,    2.0371,  1.1096,  5.233,  2.045,  3.508,  1.179,  2.533,  4.583,  4.205,  2.92,    5.849,  1.899,  0.315,  0.345,  4.806,  1.869,  2.4458,  0.833,  3.411,  1.083,    0.645,  0.636,  0.976,  4.267,  6.21,  0.68,  5.98,  1.3,  3.67,  1.81,  3.04,    1.76,  3.5,  4.68,  0.83,  3.98,  1.82,  2.78,  4.39,  3.47,  0.19,  1.33,  0.28,    0.49,  5.37,  2.4,  6.17,  6.04,  2.57,  1.71,  1.78,  0.59,  0.44,  2.74,  3.16]
cdef double* _L0C = [0.0,  6283.07585,  12566.1517,  5753.3849,  3.5231,  77713.7715,  7860.4194,  3930.2097,    11506.7698,  529.691,  1577.3435,  5884.927,  26.298,  398.149,  5223.694,  5507.553,    18849.228,  775.523,  0.067,  11790.629,  796.298,  10977.079,  5486.778,  2544.314,    5573.143,  6069.777,  213.299,  2942.463,  20.775,  0.98,  4694.003,  15720.839,  7.114,    2146.17,  155.42,  161000.69,  6275.96,  71430.7,  17260.15,  12036.46,  5088.63,    3154.69,  801.82,  9437.76,  8827.39,  7084.9,  6286.6,  14143.5,  6279.55,  12139.55,    1748.02,  5856.48,  1194.45,  8429.24,  19651.05,  10447.39,  10213.29,  1059.38,  2352.87,  6812.77,  17789.85,  83996.85,  1349.87,  4690.48]

#Earth Heliocentric Longitude coefficients (B0 and B1 in paper)
cdef int _B1L = 2
cdef double* _B1A = [9,6]
cdef double* _B1B = [3.9,1.73]
cdef double* _B1C = [5507.55,5223.69]

cdef int _B0L = 5
cdef double* _B0A = [280, 102, 80, 44, 32]
cdef double* _B0B = [3.199, 5.422, 3.88, 3.7, 4.0]
cdef double* _B0C = [84334.662, 5507.553, 5223.69, 2352.87, 1577.34]

#Earth Heliocentric Radius coefficients (R0, R1, R2, R3, R4)
cdef int _R4L = 1
cdef double* _R4A = [4]
cdef double* _R4B = [2.56]
cdef double* _R4C = [6283.08]

cdef int _R3L = 2
cdef double* _R3A = [145, 7]
cdef double* _R3B = [4.273, 3.92]
cdef double* _R3C = [6283.076, 12566.15]

cdef int _R2L = 6
cdef double* _R2A = [4359, 124, 12, 9, 6, 3]
cdef double* _R2B = [5.7846, 5.579, 3.14, 3.63, 1.87, 5.47]
cdef double* _R2C = [6283.0758, 12566.152, 0.0, 77713.77, 5573.14, 18849]

cdef int _R1L = 10
cdef double* _R1A = [103019, 1721, 702, 32, 31, 25, 18, 10, 9, 9]
cdef double* _R1B = [1.10749, 1.0644, 3.142, 1.02, 2.84, 1.32, 1.42, 5.91, 1.42, 0.27]
cdef double* _R1C = [6283.07585,  12566.1517,  0.0,  18849.23,  5507.55,  5223.69,  1577.34,  10977.08,  6275.96,  5486.78]

cdef int _R0L = 40
cdef double* _R0A = [100013989,  1670700,  13956,  3084,  1628,  1576,  925,  542,  472,  346,  329,  307,  243,  212,  186,  175,  110,  98,  86,  86,  85,  63,  57,  56,  49,  47,  45,  43,  39,  38,  37,  37,  36,  35,  33,  32,  32,  28,  28,  26]
cdef double* _R0B = [0.0,  3.0984635,  3.05525,  5.1985,  1.1739,  2.8469,  5.453,  4.564,  3.661,  0.964,  5.9,  0.299,  4.273,  5.847,  5.022,  3.012,  5.055,  0.89,  5.69,  1.27,  0.27,  0.92,  2.01,  5.24,  3.25,  2.58,  5.54,  6.01,  5.36,  2.39,  0.83,  4.9,  1.67,  1.84,  0.24,  0.18,  1.78,  1.21,  1.9,  4.59]
cdef double* _R0C = [0.0, 6283.07585, 12566.1517, 77713.7715, 5753.3849, 7860.4194, 11506.77, 3930.21, 5884.927, 5507.553, 5223.694, 5573.143, 11790.629, 1577.344, 10977.079, 18849.228, 5486.778, 6069.78, 15720.84, 161000.69, 17260.15, 529.69, 83996.85, 71430.7, 2544.31, 775.52, 9437.76, 6275.96, 4694, 8827.39, 19651.05, 12139.55, 12036.46, 2942.46, 7084.9, 5088.63, 398.15, 6286.6, 6279.55, 10447.39]

cdef double heliocentric_longitude(double jme):
    """Compute the Earth Heliocentric Longitude (L) in degrees given the Julian Ephemeris Millennium"""
    #L5, ..., L0
    cdef double L5,L4,L3,L2,L1,L0,L
    L5=L4=L3=L2=L1=L0=0
    for i in range(_L5L):
        L5 += _L5A[i]*cos(_L5B[i]+_L5C[i]*jme)
    for i in range(_L4L):
        L4 += _L4A[i]*cos(_L4B[i]+_L4C[i]*jme)
    for i in range(_L3L):
        L3 += _L3A[i]*cos(_L3B[i]+_L3C[i]*jme)
    for i in range(_L2L):
        L2 += _L2A[i]*cos(_L2B[i]+_L2C[i]*jme)
    for i in range(_L1L):
        L1 += _L1A[i]*cos(_L1B[i]+_L1C[i]*jme)
    for i in range(_L0L):
        L0 += _L0A[i]*cos(_L0B[i]+_L0C[i]*jme)
    L = polyval([L5,L4,L3,L2,L1,L0],jme)/1e8
    L = rad2deg(L) % 360
    return L
cdef double heliocentric_latitude(double jme):
    """Compute the Earth Heliocentric Latitude (B) in degrees given the Julian Ephemeris Millennium"""
    cdef double B1,B0,B
    B1=B0=0.0
    for i in range(_B1L):
        B1 += _B1A[i]*cos(_B1B[i]+_B1C[i]*jme)
    for i in range(_B0L):
        B0 += _B0A[i]*cos(_B0B[i]+_B0C[i]*jme)
    B = polyval([B1,B0], jme) / 1e8
    B = rad2deg(B) % 360
    return B
cdef double heliocentric_radius(double jme):
    """Compute the Earth Heliocentric Radius (R) in astronimical units given the Julian Ephemeris Millennium"""
    cdef double R4,R3,R2,R1,R0,R
    R4=R3=R2=R1=R0=0.0
    for i in range(_R4L):
        R4 += _R4A[i]*cos(_R4B[i]+_R4C[i]*jme)
    for i in range(_R3L):
        R3 += _R3A[i]*cos(_R3B[i]+_R3C[i]*jme)
    for i in range(_R2L):
        R2 += _R2A[i]*cos(_R2B[i]+_R2C[i]*jme)
    for i in range(_R1L):
        R1 += _R1A[i]*cos(_R1B[i]+_R1C[i]*jme)
    for i in range(_R0L):
        R0 += _R0A[i]*cos(_R0B[i]+_R0C[i]*jme)
    R = polyval([R4,R3,R2,R1,R0], jme) / 1e8
    return R
cdef void heliocentric_position(double jme, double* L_ret, double* B_ret, double* R_ret):
    """Compute the Earth Heliocentric Longitude, Latitude, and Radius given the Julian Ephemeris Millennium
        Returns (L, B, R) where L = longitude in degrees, B = latitude in degrees, and R = radius in astronimical units
    """
    L_ret[0], B_ret[0], R_ret[0] = heliocentric_longitude(jme), heliocentric_latitude(jme), heliocentric_radius(jme)

#Nutation Longitude and Obliquity coefficients (a,b,c,d,Y)
cdef int _NLOabL = 63
cdef int _NLOcdL = 49
cdef double *_NLOa = [-171996, -13187, -2274, 2062, 1426, 712, -517, -386, -301, 217, -158, 129, 123, 63, 63, -59, -58, -51, 48, 46, -38, -31, 29, 29, 26, -22, 21, 17, 16, -16, -15, -13, -12, 11, -10, -8, 7, -7, -7, -7, 6, 6, 6, -6, -6, 5, -5, -5, -5, 4, 4, 4, -4, -4, -4, 3, -3, -3, -3, -3, -3, -3, -3]
cdef double *_NLOb = [-174.2, -1.6, -0.2, 0.2, -3.4, 0.1, 1.2, -0.4, 0, -0.5, 0, 0.1, 0, 0, 0.1, 0, -0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.1, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
cdef double *_NLOc = [92025, 5736, 977, -895, 54, -7, 224, 200, 129, -95, 0, -70, -53, 0, -33, 26, 32, 27, 0, -24, 16, 13, 0, -12, 0, 0, -10, 0, -8, 7, 9, 7, 6, 0, 5, 3, -3, 0, 3, 3, 0, -3, -3, 3, 3, 0, 3, 3, 3 ]
cdef double *_NLOd = [8.9, -3.1, -0.5, 0.5, -0.1, 0, -0.6, 0, -0.1, 0.3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#N.B. Y is a 63x5 array that's been flattened
cdef double *_NLOY = [0,   0,   0,   0,   1,      -2,  0,   0,   2,   2,      0,   0,   0,   2,   2, 
                      0,   0,   0,   0,   2,      0,   1,   0,   0,   0,      0,   0,   1,   0,   0, 
                      -2,  1,   0,   2,   2,      0,   0,   0,   2,   1,      0,   0,   1,   2,   2, 
                      -2,  -1,  0,   2,   2,      -2,  0,   1,   0,   0,      -2,  0,   0,   2,   1, 
                      0,   0,   -1,  2,   2,      2,   0,   0,   0,   0,      0,   0,   1,   0,   1, 
                      2,   0,   -1,  2,   2,      0,   0,   -1,  0,   1,      0,   0,   1,   2,   1, 
                      -2,  0,   2,   0,   0,      0,   0,   -2,  2,   1,      2,   0,   0,   2,   2, 
                      0,   0,   2,   2,   2,      0,   0,   2,   0,   0,      -2,  0,   1,   2,   2, 
                      0,   0,   0,   2,   0,      -2,  0,   0,   2,   0,      0,   0,   -1,  2,   1, 
                      0,   2,   0,   0,   0,      2,   0,   -1,  0,   1,      -2,  2,   0,   2,   2, 
                      0,   1,   0,   0,   1,      -2,  0,   1,   0,   1,      0,   -1,  0,   0,   1, 
                      0,   0,   2,   -2,  0,      2,   0,   -1,  2,   1,      2,   0,   1,   2,   2, 
                      0,   1,   0,   2,   2,      -2,  1,   1,   0,   0,      0,   -1,  0,   2,   2, 
                      2,   0,   0,   2,   1,      2,   0,   1,   0,   0,      -2,  0,   2,   2,   2, 
                      -2,  0,   1,   2,   1,      2,   0,   -2,  0,   1,      2,   0,   0,   0,   1, 
                      0,   -1,  1,   0,   0,      -2,  -1,  0,   2,   1,      -2,  0,   0,   0,   1, 
                      0,   0,   2,   2,   1,      -2,  0,   2,   0,   1,      -2,  1,   0,   2,   1, 
                      0,   0,   1,   -2,  0,      -1,  0,   1,   0,   0,      -2,  1,   0,   0,   0, 
                      1,   0,   0,   0,   0,      0,   0,   1,   2,   0,      0,   0,   -2,  2,   2, 
                      -1,  -1,  1,   0,   0,      0,   1,   1,   0,   0,      0,   -1,  1,   2,   2, 
                      2,   -1,  -1,  2,   2,      0,   0,   3,   2,   2,      2,   -1,  0,   2,   2] 

cdef double ecliptic_obliquity(double jme, double delta_epsilon):
    """Calculate the true obliquity of the ecliptic (epsilon, in degrees) given the Julian Ephemeris Millennium and the obliquity"""
    cdef double u,e0,e
    u = jme/10
    e0 = polyval([2.45, 5.79, 27.87, 7.12, -39.05, -249.67, -51.38, 1999.25, -1.55, -4680.93, 84381.448], u)
    e = e0/3600.0 + delta_epsilon
    return e

cdef double dot5(double[5] x, double[5] y):
    cdef double s = 0.0
    for i in range(5):
        s += x[i]*y[i]
    return s

cdef void nutation_obliquity(double jce, double* dp_ret, double* e_ret):
    """compute the nutation in longitude (delta_psi) and the true obliquity (epsilon) given the Julian Ephemeris Century"""
    cdef double x[5]
    cdef double dp,de,e
    cdef int i
    #mean elongation of the moon from the sun, in radians:
    #x0 = 297.85036 + 445267.111480*jce - 0.0019142*(jce**2) + (jce**3)/189474
    x[0] = deg2rad(polyval([1./189474, -0.0019142, 445267.111480, 297.85036],jce))
    #mean anomaly of the sun (Earth), in radians:
    x[1] = deg2rad(polyval([-1/3e5, -0.0001603, 35999.050340, 357.52772], jce))
    #mean anomaly of the moon, in radians:
    x[2] = deg2rad(polyval([1./56250, 0.0086972, 477198.867398, 134.96298], jce))
    #moon's argument of latitude, in radians:
    x[3] = deg2rad(polyval([1./327270, -0.0036825, 483202.017538, 93.27191], jce))
    #Longitude of the ascending node of the moon's mean orbit on the ecliptic
    # measured from the mean equinox of the date, in radians
    x[4] = deg2rad(polyval([1./45e4, 0.0020708, -1934.136261, 125.04452], jce))

    dp = 0.0
    for i in range(_NLOabL):
        dp += (_NLOa[i] + _NLOb[i]*jce)*sin(dot5(x, <double*>(_NLOY+5*i)))
    dp = rad2deg(dp)/36e6

    de = 0.0
    for i in range(_NLOcdL):
        de += (_NLOc[i] + _NLOd[i]*jce)*cos(dot5(x, <double*>(_NLOY+5*i)))
    
    de = rad2deg(de)/36e6

    e = ecliptic_obliquity(julian_millennium(jce), de)
    dp_ret[0] = dp
    e_ret[0] = e

cdef double abberation_correction(double R):
    """Calculate the abberation correction (delta_tau, in degrees) given the Earth Heliocentric Radius (in AU)"""
    return -20.4898/(3600*R)

cdef void sun_longitude(double helio_L,double helio_B,double helio_R, double delta_psi, double* lambd_ret, double* beta_ret):
    """Calculate the apparent sun longitude (lambda, in degrees) and geocentric longitude (beta, in degrees) given the earth heliocentric position and delta_psi"""
    cdef double theta,beta,ll
    theta = helio_L + 180 #geocentric latitude
    beta = -helio_B
    ll = theta + delta_psi + abberation_correction(helio_R)
    lambd_ret[0] = ll
    beta_ret[0] = beta

cdef double greenwich_sidereal_time(double jd,double delta_psi,double epsilon):
    """Calculate the apparent Greenwich sidereal time (v, in degrees) given the Julian Day"""
    cdef double jc,v0,v
    jc = julian_century(jd)
    #mean sidereal time at greenwich, in degrees:
    v0 = (280.46061837 + 360.98564736629*(jd - 2451545) + 0.000387933*(jc**2) - (jc**3)/38710000) % 360
    v = v0 + delta_psi*cos(deg2rad(epsilon))
    return v

cdef void sun_ra_decl(double llambda, double epsilon, double beta, double* alpha_ret, double* delta_ret):
    """Calculate the sun's geocentric right ascension (alpha, in degrees) and declination (delta, in degrees)"""
    cdef double l,e,b,alpha,delta
    l, e, b = deg2rad(llambda), deg2rad(epsilon), deg2rad(beta)
    alpha = atan2(sin(l)*cos(e) - tan(b)*sin(e), cos(l)) #x1 / x2
    alpha = rad2deg(alpha) % 360
    delta = asin(sin(b)*cos(e) + cos(b)*sin(e)*sin(l))
    delta = rad2deg(delta)
    alpha_ret[0] = alpha
    delta_ret[0] = delta
    
cdef void sun_topo_ra_decl_hour(double latitude, double longitude, double elevation, double jd, double delta_t, double* alpha_prime_ret, double* delta_prime_ret, double* H_prime_ret):
    """Calculate the sun's topocentric right ascension (alpha'), declination (delta'), and hour angle (H')"""
    cdef double jde,jce,jme,phi,sigma,E,xi,u,x,y,delta_psi,epsilon
    cdef double llambda,beta,alpha,delta,v,H,Hr,dr,dar,delta_alpha
    cdef double alpha_prime,delta_prime,H_prime, helio_L, helio_B, helio_R
    jde = julian_ephemeris_day(jd, delta_t)
    jce = julian_century(jde)
    jme = julian_millennium(jce)

    heliocentric_position(jme,&helio_L,&helio_B,&helio_R)

    phi, sigma, E = latitude, longitude, elevation
    #equatorial horizontal parallax of the sun, in radians
    xi = deg2rad(8.794/(3600*helio_R)) #
    #rho = distance from center of earth in units of the equatorial radius
    #phi-prime = geocentric latitude
    #NB: These equations look like their based on WGS-84, but are rounded slightly
    # The WGS-84 reference ellipsoid has major axis a = 6378137 m, and flattening factor 1/f = 298.257223563
    # minor axis b = a*(1-f) = 6356752.3142 = 0.996647189335*a
    u = atan(0.99664719*tan(phi)) #
    x = cos(u) + E*cos(phi)/6378140 #rho sin(phi-prime)
    y = 0.99664719*sin(u) + E*sin(phi)/6378140 #rho cos(phi-prime)

    nutation_obliquity(jce, &delta_psi, &epsilon) #

    sun_longitude(helio_L,helio_B,helio_R, delta_psi, &llambda, &beta) #
    
    sun_ra_decl(llambda, epsilon, beta, &alpha, &delta) #

    v = greenwich_sidereal_time(jd, delta_psi, epsilon) #

    H = v + longitude - alpha #
    Hr, dr = deg2rad(H),deg2rad(delta)

    dar = atan2(-x*sin(xi)*sin(Hr), cos(dr)-x*sin(xi)*cos(Hr))
    delta_alpha = rad2deg(dar) #
    
    alpha_prime_ret[0] = alpha + delta_alpha #
    delta_prime_ret[0] = rad2deg(atan2((sin(dr) - y*sin(xi))*cos(dar), cos(dr) - y*sin(xi)*cos(Hr))) #
    H_prime_ret[0] = H - delta_alpha #

cdef void sun_topo_azimuth_zenith(double latitude, double delta_prime, double H_prime, double temperature, double pressure, double* Phi_ret, double* zenith_ret):
    """Compute the sun's topocentric azimuth and zenith angles
    azimuth is measured eastward from north, zenith from vertical
    temperature = average temperature in C (default is 14.6 = global average in 2013)
    pressure = average pressure in mBar (default 1013 = global average)
    """
    cdef double phi,dr,Hr,P,T,e0,tmp,delta_e,e,zenith,gamma,Phi
    phi = deg2rad(latitude)
    dr, Hr = deg2rad(delta_prime), deg2rad(H_prime)
    P, T = pressure, temperature
    e0 = rad2deg(asin(sin(phi)*sin(dr) + cos(phi)*cos(dr)*cos(Hr)))
    tmp = deg2rad(e0 + 10.3/(e0+5.11))
    delta_e = (P/1010.0)*(283.0/(273+T))*(1.02/(60*tan(tmp)))
    e = e0 + delta_e
    zenith_ret[0] = 90 - e

    gamma = rad2deg(atan2(sin(Hr), cos(Hr)*sin(phi) - tan(dr)*cos(phi))) % 360
    Phi_ret[0] = (gamma + 180) % 360 #azimuth from north

cdef void norm_lat_lon(double* lat, double* lon):
    A,O = lat[0],lon[0]
    if A < -90 or A > 90:
        #convert to cartesian and back
        x = cos(deg2rad(O))*cos(deg2rad(A))
        y = sin(deg2rad(O))*cos(deg2rad(A))
        z = sin(deg2rad(A))
        r = sqrt(x**2 + y**2 + z**2)
        O = rad2deg(atan2(y,x)) % 360
        A = rad2deg(asin(z/r))
    elif O < 0 or O > 360:
        O = O % 360
    lat[0],lon[0] = A,O

cdef void pos(t,double lat,double lon,double elev,double temp,double press,double dt, double* ret):
    """Compute azimute,zenith,RA,dec,H all in degrees"""
    cdef double jd,RA,dec,H,azimuth,zenith
    norm_lat_lon(&lat,&lon)
    jd = cjulian_day(t)
    sun_topo_ra_decl_hour(lat, lon, elev, jd, dt, &RA, &dec, &H)
    sun_topo_azimuth_zenith(lat, dec, H, temp, press, &azimuth, &zenith)
    ret[0] = azimuth
    ret[1] = zenith

cdef void topo_pos(t,double lat,double lon,double elev,double temp,double press,double dt, double* ret):
    """compute RA,dec,H, all in degrees"""
    cdef double jd,RA,dec,H
    norm_lat_lon(&lat,&lon)
    jd = cjulian_day(t)
    sun_topo_ra_decl_hour(lat, lon, elev, jd, dt, &RA, &dec, &H)
    ret[0] = RA
    ret[1] = dec
    ret[2] = H

cdef void full_pos(t,double lat,double lon,double elev,double temp,double press,double dt, double* ret):
    """Compute azimute,zenith,RA,dec,H all in degrees"""
    cdef double jd,RA,dec,H,azimuth,zenith
    norm_lat_lon(&lat,&lon)
    jd = cjulian_day(t)
    sun_topo_ra_decl_hour(lat, lon, elev, jd, dt, &RA, &dec, &H)
    sun_topo_azimuth_zenith(lat, dec, H, temp, press, &azimuth, &zenith)
    ret[0] = azimuth
    ret[1] = zenith
    ret[2] = RA
    ret[3] = dec
    ret[4] = H

def julian_day(dt):
    """julian_day(dt)
    Convert UTC datetimes or UTC timestamps to Julian days

    Parameters
    ----------
    dt : array_like
        UTC datetime objects or UTC timestamps (as per datetime.utcfromtimestamp)

    Returns
    -------
    jd : ndarray
        datetimes converted to fractional Julian days
    """
    dts = np.array(dt)
    if len(dts.shape) == 0:
        return cjulian_day(dt)

    jds = np.empty(dts.shape)
    for i,d in enumerate(dts.flat):
        jds.flat[i] = cjulian_day(d)
    return jds

def arcdist(p0,p1,radians=False):
    """arcdist(p0,p1,radians=False)
    Angular distance between azimuth,zenith pairs
    
    Parameters
    ----------
    p0 : array_like, shape (..., 2)
    p1 : array_like, shape (..., 2)
        p[...,0] = azimuth angles, p[...,1] = zenith angles
    radians : boolean (default False)
        If False, angles are in degrees, otherwise in radians

    Returns
    -------
    ad :  array_like, shape is broadcast(p0,p1).shape
        Arcdistances between corresponding pairs in p0,p1
        In degrees by default, in radians if radians=True
    """
    #formula comes from translating points into cartesian coordinates
    #taking the dot product to get the cosine between the two vectors
    #then arccos to return to angle, and simplify everything assuming real inputs
    p0,p1 = np.array(p0), np.array(p1)
    a0,z0 = p0[...,0], p0[...,1]
    a1,z1 = p1[...,0], p1[...,1]
    d = np.empty(np.broadcast(a0,a1).shape)
    cdef np.broadcast it = np.broadcast(a0,z0,a1,z1,d)
    cdef bint r = radians
    cdef double _a0,_z0,_a1,_z1,_d
    while np.PyArray_MultiIter_NOTDONE(it):
        _a0 = (<double*>np.PyArray_MultiIter_DATA(it,0))[0]
        _z0 = (<double*>np.PyArray_MultiIter_DATA(it,1))[0]
        _a1 = (<double*>np.PyArray_MultiIter_DATA(it,2))[0]
        _z1 = (<double*>np.PyArray_MultiIter_DATA(it,3))[0]
        if not r:
            _a0,_z0 = deg2rad(_a0),deg2rad(_z0)
            _a1,_z1 = deg2rad(_a1),deg2rad(_z1)
        _d = acos(cos(_z0)*cos(_z1)+cos(_a0-_a1)*sin(_z0)*sin(_z1))
        if not r:
            _d = rad2deg(_d)
        (<double*>np.PyArray_MultiIter_DATA(it,4))[0] = _d
        np.PyArray_MultiIter_NEXT(it)
    return d

def sunpos_az(dt, latitude, longitude, elevation, temperature=None, pressure=None, delta_t=0, radians=False):
    """sunpos_az(dt,latitude,longitude,elevation,temperature=None,pressure=None,delta_t=0,radians=False)
    Compute the azimuth and zenith angles of the sun as viewed at the given time and location.

    Parameters
    ----------
    dt : array_like
        UTC datetime objects or UTC timestamps (as per datetime.utcfromtimestamp) representing the times of observations
    latitude, longitude : array_like
        decimal degrees, positive for north of the equator and east of Greenwich
    elevation : array_like
        meters, relative to the WGS-84 ellipsoid
    temperature : array_like or None, optional
        celcius, default is 14.6 (global average in 2013)
    pressure : array_like or None, optional
        millibar, default is 1013 (global average in ??)
    delta_t : array_like, optional
        seconds, default is 0, difference between the earth's rotation time (TT) and universal time (UT)
    radians : {True, False}, optional
        return results in radians if True, degrees if False (default)

    Returns
    -------
    coords : ndarray, (...,2)
        The shape of the array is parameters broadcast together, plus a final dimension for the coordinates.
        coords[...,0] = observed azimuth angle, measured eastward from north
        coords[...,1] = observed zenith angle, measured down from vertical
    """
    if temperature is None:
        temperature = 14.6
    if pressure is None:
        pressure = 1013
    
    #6367444 = radius of earth
    #numpy broadcasting
    b = np.broadcast(dt,latitude,longitude,elevation,temperature,pressure,delta_t)
    res = np.empty((b.size,2))
    cdef np.double_t[:,:] res_view = res
    cdef size_t i
    for i,x in enumerate(b):
        t,lat,lon,el,temp,pres,delt = x
        pos(t,lat,lon,el,temp,pres,delt,&res_view[i,0])
    if radians:
        res = np.deg2rad(res)
    return res.reshape(b.shape+(2,))

def sunpos_adh(dt, latitude, longitude, elevation, temperature=None, pressure=None, delta_t=0, radians=False):
    """sunpos_adh(dt, latitude, longitude, elevation, temperature=None, pressure=None, delta_t=0, radians=False)
    Compute the right ascension, declination, and hour angles of the sun as viewed at the given time and location.

    Parameters
    ----------
    dt : array_like
        UTC datetime objects or UTC timestamps (as per datetime.utcfromtimestamp) representing the times of observations
    latitude, longitude : array_like
        decimal degrees, positive for north of the equator and east of Greenwich
    elevation : array_like
        meters, relative to the WGS-84 ellipsoid
    temperature : array_like or None, optional
        celcius, default is 14.6 (global average in 2013)
    pressure : array_like or None, optional
        millibar, default is 1013 (global average in ??)
    delta_t : array_like, optional
        seconds, default is 0, difference between the earth's rotation time (TT) and universal time (UT)
    radians : {True, False}, optional
        return results in radians if True, degrees if False (default)

    Returns
    -------
    coords : ndarray, (...,3)
        The shape of the array is parameters broadcast together, plus a final dimension for the coordinates.
        coords[...,0] = topocentric right ascension
        coords[...,1] = topocentric declination
        coords[...,2] = topocentric hour angle
    """
    if temperature is None:
        temperature = 14.6
    if pressure is None:
        pressure = 1013
    
    #6367444 = radius of earth
    #numpy broadcasting
    b = np.broadcast(dt,latitude,longitude,elevation,temperature,pressure,delta_t)
    res = np.empty((b.size,3))
    cdef np.double_t[:,:] res_view = res
    cdef size_t i
    for i,x in enumerate(b):
        t,lat,lon,el,temp,pres,delt = x
        topo_pos(t,lat,lon,el,temp,pres,delt,&res_view[i,0])
    if radians:
        res = np.deg2rad(res)
    return res.reshape(b.shape+(3,))

def sunpos(dt, latitude, longitude, elevation, temperature=None, pressure=None, delta_t=0, radians=False):
    """sunpos(dt, latitude, longitude, elevation, temperature=None, pressure=None, delta_t=0, radians=False)
    Compute the observed and topocentric coordinates of the sun as viewed at the given time and location.

    Parameters
    ----------
    dt : array_like
        UTC datetime objects or UTC timestamps (as per datetime.utcfromtimestamp) representing the times of observations
    latitude, longitude : array_like
        decimal degrees, positive for north of the equator and east of Greenwich
    elevation : array_like
        meters, relative to the WGS-84 ellipsoid
    temperature : array_like or None, optional
        celcius, default is 14.6 (global average in 2013)
    pressure : array_like or None, optional
        millibar, default is 1013 (global average in ??)
    delta_t : array_like, optional
        seconds, default is 0, difference between the earth's rotation time (TT) and universal time (UT)
    radians : {True, False}, optional
        return results in radians if True, degrees if False (default)

    Returns
    -------
    coords : ndarray, (...,5)
        The shape of the array is parameters broadcast together, plus a final dimension for the coordinates.
        coords[...,0] = observed azimuth angle, measured eastward from north
        coords[...,1] = observed zenith angle, measured down from vertical
        coords[...,2] = topocentric right ascension
        coords[...,3] = topocentric declination
        coords[...,4] = topocentric hour angle
    """

    if temperature is None:
        temperature = 14.6
    if pressure is None:
        pressure = 1013
    
    #6367444 = radius of earth
    #numpy broadcasting
    b = np.broadcast(dt,latitude,longitude,elevation,temperature,pressure,delta_t)
    res = np.empty((b.size,5))
    cdef np.double_t[:,:] res_view = res
    cdef size_t i
    for i,x in enumerate(b):
        t,lat,lon,el,temp,pres,delt = x
        full_pos(t,lat,lon,el,temp,pres,delt,&res_view[i,0])
    if radians:
        res = np.deg2rad(res)
    return res.reshape(b.shape+(5,))

def citation():
    return 'I. Reda and A. Andreas, "Solar position algorithm for solar radiation applications," Solar Energy, vol. 76, no. 5, pp. 577-589, 2004'

def bibtex():
    return '@article{RN120,\n    author = {Reda, Ibrahim and Andreas, Afshin},\n    title = {Solar position algorithm for solar radiation applications},\n    journal = {Solar Energy},\n    volume = {76},\n    number = {5},\n    pages = {577-589},\n    ISSN = {0038-092X},\n    DOI = {10.1016/j.solener.2003.12.003},\n    url = {http://www.sciencedirect.com/science/article/pii/S0038092X0300450X},\n    year = {2004},\n    type = {Journal Article}\n}'

