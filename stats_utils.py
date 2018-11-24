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
from scipy import stats

def lanczos(a,n):
    """lanczos filter kernel, order a, n points"""
    x = linspace(-a+.5,a-.5,n)
    y = sinc(x)*sinc(x/a)
    return y/sum(y)

def bin(x,y,bins,xmin=None,xmax=None,min_count=1,keep_oob=True,keep_empty=True):
    """x and y are 1d ndarrays"""
    x = asarray(x)
    y = asarray(y)
    if isscalar(bins):
        #a number of bins, not the bins 
        if xmin is None:
            xmin = x.min()
        if xmax is None:
            xmax = x.max()
        nbins = bins
        #set up bins for digitize
        w = (xmax-xmin)/nbins
        bins = linspace(xmin,xmax,nbins+1) #includes endpoint
        #bin midpoints
        mids = bins[:-1] + (w/2)
    else:
        nbins = len(bins)-1
        mids = (bins[:-1] + bins[1:])/2
    #digitize x -- get bin index for each x
    xidx = digitize(x,bins)
    #send ys to the appropriate bin
    y_binned = []
    not_empty = []
    if keep_oob:
        rng = range(nbins+2) #include out-of-bounds bins
    else:
        rng = range(1,nbins+1) #don't include out-of-bounds bins
    
    for i in rng:
        y_b = y[xidx==i]
        if len(y_b) >= min_count:
            y_binned.append(y_b)
            if not keep_empty:
                #we're not keeping empty bins, so the non-empty mids need to be tracked
                if i == 0: #the 0th bin is all numbers less than xmin
                    not_empty.append(-float('inf'))
                elif i == nbins+1: #the last bin is all greater than xmax
                    not_empty.append(float('inf'))
                else: 
                    not_empty.append(mids[i-1])
        elif keep_empty:
            y_binned.append(array([]))
    if keep_empty:
        return mids, y_binned
    else:
        return array(not_empty), y_binned

def wrap(x,vmin=0.0,vmax=1.0):
    return (x-vmin)%(vmax-vmin)+vmin

def clamp(x,vmin=0.0, vmax=1.0):
    if isscalar(x):
        if vmax is not None:
            x = min(x,vmax)
        if vmin is not None:
            x = max(x,vmin)
        return x
    else:
        x = asarray(x)
        x_c = copy(x)
        if vmax is not None:
            x_c[x >= vmax] = vmax
        if vmin is not None:
            x_c[x <= vmin] = vmin
        return x_c

def bin_wrap(x,y,bins,xmin=None,xmax=None,center=0,min_count=1,keep_empty=True):
    #center is a relative offset of the 1st bin's center from xmin
    #center should be in [0,0.5]
    #center = 0 -> 1st bin is centered on xmin
    #center = 0.5 -> 1st bin's left edge is at xmin
    x = asarray(x)
    y = asarray(y)
    if isscalar(bins):
        nbins = bins
        if xmin is None:
            xmin = x.min()
        if xmax is None:
            xmax = x.max()
        w = (xmax-xmin)/nbins
        bins = linspace(xmin,xmax,nbins+1) + w*(center-0.5)
    
    bmin = bins[0]
    bmax = bins[-1]
    brange = bmax - bmin
    x_wrapped = ((x-bmin) % brange)+bmin
    mids,y_binned = bin(x_wrapped,y,bins,min_count=min_count,keep_oob=False,keep_empty=keep_empty)
    #because the data was wrapped, we don't need the out-of-bounds bins
    return mids,y_binned

def likelihood_same_mean(mahal_dist,ndims=1):
    """What is the likelihood that you drew x from N(mu,sigma)
    mahal_dist = Mahalanobis distance (x-mu)/sigma or sqrt( (x-mu).T.dot(inv(S).dot((x-mu))), S is covariance matrix
    ndims = number of dimensions of x
    """
    return stats.chi2.cdf(mahal_dist**2,ndims)

def mahal_distance_for_likelihood(likelihood,ndims=1):
    return sqrt(stats.chi2.isf(1-likelihood,ndims))

def rms(x,axis=None,keepdims=False):
    return sqrt(mean(x**2,axis,keepdims=keepdims))

def mean_conf(data,conf=0.90):
    n = len(data)
    #mean, and std error of mean
    m,se = np.mean(data), stats.sem(data)
    #get the interval from Student's t dist
    iv = stats.t.interval(conf,n-1,scale=se)
    return m,iv

def angle_diff(a1,a2,period=2*pi):
    """(a1 - a2 + d) % (2*d) - d; d = period/2"""
    d = period/2
    return ((a1 - a2 + d) % (period)) - d

def angle_mean(ang,axis=None,period=2*pi):
    """returns the circular mean of angles"""
    #uses the 1st angular moment:
    a=2*pi/period
    m1 = np.mean(np.exp(1j*ang*a),axis=axis)
    return np.angle(m1)/a

def angle_std(ang,axis=None,period=2*pi):
    """Returns the circular standard deviation of angles"""
    a = 2*pi/period
    m1 = np.mean(np.exp(1j*ang*a),axis=axis)
    return np.sqrt(-2*np.log(np.abs(m1)))/a

def angle_var(ang,axis=None,period=2*pi):
    a = 2*pi/period
    m1 = np.mean(np.exp(1j*ang*a),axis=axis)
    return -2*np.log(np.abs(m1))/a**2

def angle_mean_std(ang,axis=None,period=2*pi):
    """returns the circular mean and standard deviation of angles"""
    #take the 1st angular moment
    # nth moment: m_n = mean(exp(1j*ang)**n)
    a = 2*pi/period
    m1 = np.mean(np.exp(1j*ang*a),axis=axis)
    mean = np.angle(m1)/a
    std = np.sqrt(-2*np.log(np.abs(m1)))/a
    return mean, std

def angle_rms(ang,axis=None,period=2*pi):
    """returns the rms of angles, uses the property that rms(x)**2 = mean(x)**2 + std(x)**2"""
    #rms(x)**2 = mean(x)**2 + std(x)**2
    #sqrt(E[X**2]) = E[X]**2 + sqrt(E[(X - E[X])**2])
    m,s = angle_mean_std(ang,axis,period)
    return hypot(m, s)

def angle_wrap(p,period=2*pi):
    if isscalar(period):
        vmin = -period/2
        vmax = period/2
    else:
        vmin = period[0]
        vmax = period[1]
    return (p-vmin)%(vmax-vmin)+vmin

def angle_unwrap(p,axis=-1,period=2*pi):
    """angle_unwrap(p,axis=-1,period=2*pi)
    Unwrap angles by changing deltas between values to their minimum relative to the period
    
    Parameters
    ----------
    p : array_like
        Angles to be unwrapped
    period : float, default = 2*pi
        The maximum angle. Differences between subsequent p are set to their minimum relative to this value. Set to 2*pi for radians, 360 for degrees, 1 for turns, etc.
    axis : int, default = -1
        Axis along which unwrap will be performed.
    
    Returns
    -------
    out : ndarray
        unwrapped angles
    """
    p = asarray(p)
    nd = len(p.shape)
    dd = diff(p, axis=axis)
    slice1 = [slice(None, None)]*nd     # full slices
    slice1[axis] = slice(1, None)
    hc = period/2
    ddmod = mod(dd + hc, period) - hc
    copyto(ddmod, hc, where=(ddmod == -hc) & (dd > 0))
    ph_correct = ddmod - dd
    copyto(ph_correct, 0, where=abs(dd) < hc)
    up = array(p, copy=True, dtype='d')
    up[slice1] = p[slice1] + ph_correct.cumsum(axis)
    return up

def circle_mask(radius,size=None,offset=None,inner=0,subsample_limit=4,center=False):
    def subsample(x,y,sz,r,lim):
        d = hypot(x,y)
        if lim==0: #hit recursion limit
            #return area if x,y is inside circle
            return sz**2 if d < r else 0.0
        elif d + 0.70711*sz < r: #totally inside circle
            return sz**2
        elif d - 0.70711*sz > r: #totally outside circle
            return 0.0
        else: #on edge, recurse into quadrants
            s,o = sz/2, sz/4
            return subsample(x+o,y+o,s,r,lim-1) + \
                   subsample(x+o,y-o,s,r,lim-1) + \
                   subsample(x-o,y-o,s,r,lim-1) + \
                   subsample(x-o,y+o,s,r,lim-1)
    if offset is None:
        y0,x0 = 0,0
    else:
        y0,x0 = offset
    if size is None:
        size=2*radius+1
    if isscalar(size):
        size = (size,size)
    if center:
        y0 += 0.5*size[0]-0.5-radius
        x0 += 0.5*size[1]-0.5-radius
    coeffs = empty(size)
    for r in range(size[0]):
        for c in range(size[1]):
            x,y = c-radius,r-radius
            coeffs[r,c] = subsample(x-x0,y-y0,1,radius,subsample_limit)
            if inner > 0:   
                coeffs[r,c] -= subsample(x-x0,y-y0,1,inner,subsample_limit) 
    return coeffs
