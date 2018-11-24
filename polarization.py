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
import numpy as np

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

class Mueller:
    @staticmethod
    def rotator(angle):
        """Rotates Stokes Vector by angle radians from S1 to S2. This is the opposite of the Polarized Light book convention"""
        c = np.cos(2*angle)
        s = np.sin(2*angle)
        return np.array([[1,0,0,0],[0,c,-s,0],[0,s,c,0],[0,0,0,1]])

    @staticmethod
    def rotate(mat,angle):
        """Rotate a Mueller matrix from S1 to S2 by angle radians"""
        return np.dot(Mueller.rotator(angle), np.dot(mat, Mueller.rotator(-angle)))

    @staticmethod
    def polarizer(px,py,angle=0):
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
        c = np.cos(phase)
        s = np.sin(phase)
        m = np.array([[1,0,0,0],[0,1,0,0],[0,0,c,s],[0,0,-s,c]])
        if angle != 0:
            return Mueller.rotate(m,angle)
        else:
            return m

    @staticmethod
    def depolarizer(dp):
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

class Stokes:
    @staticmethod
    def intensity(s,axis=0):
        return np.take(s,0,axis)
    @staticmethod
    def aop(s,axis=0):
        s1 = np.take(s,1,axis)
        s2 = np.take(s,2,axis)
        return 0.5*np.arctan2(s2,s1)
    @staticmethod
    def dop(s,axis=0):
        s0 = np.take(s,0,axis)
        s123 = np.take(s,(1,2,3),axis)
        return np.linalg.norm(s123,axis=axis)/s0
    @staticmethod
    def dolp(s,axis=0):
        s0 = np.take(s,0,axis)
        s12 = np.take(s,(1,2),axis)
        return np.linalg.norm(s12,axis=axis)/s0
    @staticmethod
    def docp(s,axis=0):
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
        return ki-2*n*Fresnel._vdot(n,ki)

    @staticmethod
    def reflect(ki,n):
        """Reflect a plane wave with propagation vector ki off an interface with unit normal vector n.
        The normal vector points towards the incident side of the surface.
        """
        ki,n = map(np.asarray, (ki,n))
        return Fresnel._reflect(ki,n)

    @staticmethod
    def _refract(ki,n,mi,mt):
        #use inner() instead of dot() for broadcasting
        kll = ki-n*Fresnel._vdot(ki,n)
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
        tt = 2*Fresnel._vdot(ki,n,keepdims=False)
        ts = tt / Fresnel._vdot(ki + kt, n,keepdims=False)
        tp = mi*mt*tt / Fresnel._vdot(mt**2*ki + mi**2*kt, n,keepdims=False)
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
