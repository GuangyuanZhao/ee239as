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
import h5py
from pylab import *
import pdb
#model: each pixel is modeled by:
#   I_pix = dot(A,S_in) + d
#where A is the analysis vector, S_in is the incident stokes vector, and d is the offset
#note that A encapsulates the conversion gain of the pixel

#we're going to learn A and d for each pixel by using a linear regression
#这个是怎么作用？
## Functions ##
###############

def dolp(S0,S1,S2):
    return hypot(S1,S2)/S0

def aop(S1,S2):
    return 0.5*arctan2(S2,S1)

def mrot(th):
    """Rotate stokes vector by th from s1 to s2 (negative of "Polarized Light" book standard)"""
    return array([[1, 0,          0,         0],
                  [0, cos(2*th), -sin(2*th), 0],
                  [0, sin(2*th),  cos(2*th), 0],
                  [0, 0,          0,         1]])

def pa_to_s(stim):
    """Compute stokes vector given (power, angle_radians)"""
    p,a = stim
    return dot(mrot(a),[p,p,0,0])

def ap_to_s(stim):
    """Compute stokes vector given (angle_radians, power)"""
    a, p = stim
    return dot(mrot(a),[p,p,0,0])

def a_to_s(stim):
    """compute stokes vector given just (angle_radians)"""
    return dot(mrot(stim),[1,1,0,0])

def p_to_s(stim):
    return array([stim,0,0,0])

def pixel_parameters(stimulus,data,dark=True,photometric=False,stimulus_to_s=None):
    """Compute pixel parameters (A, d) for each pixel in data
    Parameters:
     stimulus has shape (nsamples, ...)
     data has shape (nsamples, ...)
     dark (True)
        set to False to not compute d
     photometric (False)
        will keep A in the units of the stokes vectors rather than normalizing to the mean of the data
     stimulus_to_s (None)
        is a function that takes a stimulus value and returns the corresponding stokes vector with shape (nstokes,)
        if None, the value of the stimulus is used directly; so stimulus has shape (nsamples, nstokes)
    Returns:
     A if dark is False
     (A,d) tuple if dark is True
        A has shape data.shape[1:] + (nstokes,)
        d has shape data.shape[1:]
    """
    if len(stimulus.shape) == 1:
        stim = stimulus[:,None]
    else:
        stim = stimulus

    nsamples = stim.shape[0]

    if stimulus_to_s is None:
        nstokes = stim.shape[1]
    else:
        nstokes = stimulus_to_s(stim[0]).shape[0]

    if dark:
        nstokes += 1 #augment for dark values

    #we're going to solve for [A,d]:
    # data = dot([S, 1],[A, d])
    #   data       is (M,K) = (nsamples) x (npix)
    #   [S,1]      is (M,N) = (nsamples) x (nstokes + 1)
    #   [A,d] will be (N,K) = (nstokes + 1) x (npix)

    #vectorize data
    data_vec = data.reshape((nsamples, -1))

    #make (possibly augmented) stokes vectors for each sample
    stokes = empty((nsamples,nstokes))

    #make a view into the unaugmented portion if necessary
    if dark:
        stokes_stim = stokes[...,:-1]
    else:
        stokes_stim = stokes

    #fill with stimulus values
    if stimulus_to_s is None:
        stokes_stim[:] = stim
    else:
        for i in range(nsamples):
            stokes_stim[i] = stimulus_to_s(stim[i])
    
    if not photometric and not normalize:
        #convert stokes vectors to data's units
        #multiply the stokes vectors by (mean of data)/(mean of s[:,0])
        stokes *= 2*mean(data)/mean(stokes[:,0])
    
    #the augmented part, if it exists, is all 1's, regardless of the photometric business
    if dark:
        stokes[:,-1] = 1
    
    #allocate & vectorize output space
    Ad = empty(data.shape[1:]+(nstokes,))
    Ad_vec = Ad.reshape((-1,nstokes)).T #transpose for lstsq
    
    #solve
    Ad_vec[:] = lstsq(stokes,data_vec)[0]
    if dark:
        A = Ad[...,:-1]
        d = Ad[...,-1]
        return (A,d)
    else:
        #not augmented, we just have A
        return Ad


def pixel_pattern(shape,A,angles=None):
    """Determine the (shape) pattern of angles in A
      shape is the shape of the pattern. len(shape) <= len(A.shape)
      if angles is given, the returned pattern will be rounded to the closest angles
      A has shape (...,nstokes)
    """
    #find the mean A[...,1] and A[...,2] for each starting point in the shape
    A_mean = empty(shape+(2,))
    for i in range(prod(shape)):
        start = unravel_index(i,shape)
        slices = tuple(slice(a,None,b) for a,b in zip(start,shape))
        A_mean[slices] = mean(A[slices][...,1:3].reshape((-1,2)),axis=0)
    #no angles to fit to: just return avg. aop of filters
    if angles is None:
        return aop(A_mean[0],A_mean[1])
    #try to find best fit angles
    #normalize
    A_mean /= norm(A_mean,axis=-1)[...,None] #the norm loses an axis so we add it back
    pattern = empty(shape)*nan
    for a in angles:
        A_desired = [cos(2*a),sin(2*a)]
        i = argmin(norm(A_mean-A_desired,axis=-1))
        pattern.flat[i] = a
    return pattern


def compute_calibration(A, slices, A_ideal, no_gain=False):
    """Compute calibration matrices C, such that dot(C,A_sliced) = A_ideal
    Parameters:
      A : shape (..., nstokes)
      slices : shape (nslices) list of slice objects
        Each slice in slices is applied to A and concatenated to create A_sliced
        For each slice s in slices, A[s] must have the same shape (slice_shape...)
      A_ideal : shape (nstokes, nslices)
      no_gain : default False
        if True, normalize each calibration matrix to have no net gain/loss
        (removes any flat-field corrections!)
    Returns:
      C : shape (slice_shape..., nslices, nslices)
    """
    #reorder A according to the slices
    #A[s][...,None] has shape (slice_shape..., nstokes, 1)
    A_sliced = concatenate([A[s][...,None] for s in slices],axis=-1)
    #A_sliced now has shape (slice_shape..., nstokes, nslices)
    slice_shape = A_sliced.shape[:-2]
    A_sliced = A_sliced.reshape((-1,) + A_sliced.shape[-2:])
    #A_sliced now has shape (nsuperpix, nstokes, nslices)

    #C's shape is (slice_shape..., nslices, nslices)
    C = empty(slice_shape+(A_sliced.shape[-1],)*2)
    C_vec = C.reshape((-1,)+C.shape[-2:])

    #the calibration eq'n is:
    # Since I_pix = dot(S, A) + d  (for I_pix, S as row vectors)
    # then S = (I_pix - d) / A (for a group of pixels w/ the same S!)
    # I_cal = dot(S, A_ideal)
    # so I_cal = dot((I_pix - d) / A, A_ideal) = dot(I_pix - d, A \ A_ideal)
    # let C = A \ A_ideal ==> I_cal = dot(I_pix - d, C) and A*C = A_ideal
    # we like to do column vectors, so store C.T so that I_cal = dot(C, I_pix - d)
    #do the regression for each superpixel:
    for ix in range(A_sliced.shape[0]):
        C_vec[ix].T[:] = linalg.lstsq(A_sliced[ix], A_ideal)[0]

    return C


## Load Data ##
###############
if __name__ == '__main__':
    print('Loading camera data...',flush=True)
    datafile = h5py.File('cam.h5','r')
    angles = deg2rad(array(datafile['angles'])) # (nangles,)
    powers = array(datafile['powers']) # (npowers,)
    data = array(datafile['data']) # (npowers,nangles,rows,cols)

    rcoffset = array([0,0],dtype=int32)

    datafile.close()

    ## Compute ##
    #############

    print('Solving for pixel parameters...',flush=True)

    #the stimulus is (angle,power) for each sample in this case
    #we have to reshape the data and make the array of angle,power pairs
    data = data.reshape((-1,)+data.shape[2:])
    a,p = meshgrid(angles,powers)
    ang_pow = array(list(zip(a.flat,p.flat)))
    A,d = pixel_parameters(ang_pow,data,True,False,ap_to_s)

    print('Mean A0, A1, A2: {}, {}, {}'.format(mean(A[...,0]),mean(A[...,1]),mean(A[...,2])))
    print('Mean d: {}'.format(mean(d)))

    print('Determining pixel pattern...',flush=True)
    order = [0,90,45,135]
    pattern = rad2deg(pixel_pattern((2,2),A,deg2rad(order)))

    print(pattern)

    print('Computing calibration matrices...',flush=True)
    #if we wanted to reconstruct data_vec from S_in and A_vec:
    #data_r = empty_like(data,dtype=float32)
    #data_r_vec = data_r.reshape(data_vec.shape)
    #data_r_vec[:] = dot(S_in, A_vec)


    starts = {pattern[p]:p for p in [(0,0),(0,1),(1,0),(1,1)]}
    #slices for vectorizing super pixels in A: first 2 dimensions are image dimensions, last dimension is individual pixel's A vector
    A_slices = [(slice(starts[a][0],None,2), slice(starts[a][1],None,2), slice(None)) for a in order]

    #what's the ideal A for those 4 pixels?
    orad = deg2rad(array(order))
    A_ideal = 0.5*array([ones_like(orad), cos(2*orad), sin(2*orad), zeros_like(orad)])

    C = compute_calibration(A, A_slices, A_ideal)


    #print('Computing reduction matrices...')

    #R = empty(Ag.shape)
    #R_vec = R.reshape((-1,)+R.shape[-2:])
    #for ix in range(Ag_vec.shape[0]):
    #    R_vec[ix].T[:] = linalg.pinv(Ag_vec[ix])

    print('writing results to disk...')
    calfile = h5py.File('cal.h5','w')
    opts = {'chunks':True, 'compression':9, 'shuffle':True, 'fletcher32':True} #dataset options

    grp = calfile.create_group('parameters')
    h5pat = grp.create_dataset('pattern',data=pattern)
    h5rc0 = grp.create_dataset('rc0',data=rcoffset)
    h5drk = grp.create_dataset('d',data=d,dtype=float32,**opts)
    grp.create_dataset('A',data=A,dtype=float32,**opts)

    grp = calfile.create_group('calibration')
    grp['darks'] = h5drk
    grp.create_dataset('gains',data=C,dtype=float32,**opts)
    grp['attributes/pattern'] = h5pat
    h5y0 = grp.create_dataset('attributes/y0',data=rcoffset[0])
    h5x0 = grp.create_dataset('attributes/x0',data=rcoffset[1])

    calfile.close()
