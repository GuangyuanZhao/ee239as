# preprocess.py
# 
# Copyright (c) 2016 Samuel B. Powell, Washington University in St Louis
# MIT License
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

from pylab import *
from scipy import signal
import h5py, os, pdb
from glob import glob
from collections import namedtuple
import numexpr as ne

def file_newer(path0,path1):
    """return true is path0 was modified after path1"""
    return os.stat(path0).st_mtime > os.stat(path1).st_mtime

def make_slices(pattern):
    # Here it is obvious that the author want to separate the images.
    #make slice objects that split the image into the different filtered components
    #ie. i0 = data[slices[0]], i90 = data[slices[1]], i45 = data[slices[2]], i135 = data[slices[3]]
    positions = [(0,0),(0,1),(1,0),(1,1)]
    starts = {pattern[p]:p for p in positions} #what angle is each start position?
    angles = [0,90,45,135]
    slices = [tuple([Ellipsis, slice(starts[a][0],None,2), slice(starts[a][1],None,2)]) for a in angles]
    return slices

CalData = namedtuple('CalData',['gains','darks','pattern','slices'])

def load_cal(fname):
    """returns gains, darks, pattern, slices"""
    if fname is not None and fname is not "":
        with h5py.File(fname,'r') as calfile:
            gains = np.array(calfile['calibration/gains'], copy=True)
            darks = np.array(calfile['calibration/darks'], copy=True)
            pattern = np.array(calfile['calibration/attributes/pattern'],copy=True)
    else:
        gains = None
        darks = None
        pattern = np.array([[0,45],[135,90]]) #default pixel layout
    slices = make_slices(pattern)
    return CalData(gains,darks,pattern,slices)

def interpolate(img_stack, slices, method='bcspline'):
    """interpolate a single frame. method can be none, bilinear, lanczos3, or bcspline"""
    
    #TODO: make this more efficient!!
    if method == 'none':
        return img_stack
    
    #create the filter based on the type of interpolation:
    if method == 'bcspline':
        #reconstruction/interpolation filter
        #3rd order bspline function
        filt = signal.bspline([-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5],3)
    elif method == 'bilinear':
        filt = [0.5,1.0,0.5]
    elif method == 'lanczos3':
        #x = linspace(-2.5,2.5,11) #evenly spaced from -2.5 to 2.5 by 0.5 increments
        #filt = sinc(x/3)*sinc(x)
        filt = [2.44565217e-02, 0, -1.35869565e-01, 0, 6.11413043e-01, 1,
                6.11413043e-01, 0, -1.35869565e-01, 0, 2.44565217e-02]

    #allocate the output array:
    #img_stack.shape is (w/2, h/2, 4)
    #we're going to expand to (w,h,4):
    s = list(img_stack.shape)
    s[:2] = np.multiply(s[:2], 2) #element-wise multiply (w,h,4)    
    #interpolated result:
    interp_stack = zeros(s)

    #interpolation coefficients for intermediate steps
    c_jk = zeros(s[:2]) #(w,h)    
    
    #loop over last dimension (i0,i90,...)
    for j in xrange(s[-1]):
        #zero the coefficients
        c_jk[...] = 0.0
        #get coefficients
        if method in ['bilinear', 'lanczos3']:
            #coefficients are just the image slices
            c_jk[slices[j]] = img_stack[...,j]
        elif method == 'bcspline':
            #coefficients are bicubic spline coefficients for slice
            c_jk[slices[j]] = signal.cspline2d(img_stack[...,j])
        #convolve (filters are seperable, so we can use sepfir2d)
        interp_stack[...,j] = signal.sepfir2d(c_jk,filt,filt)
    #return interpolated image:
    return interp_stack

def process(data, caldata, imethod='none'):
    ## Calibration ##
    gains,darks,pattern,slices = caldata
    #first subtract dark values
    if darks is not None:
        d = ne.evaluate('data-darks') #broadcasts if data is (...,r,c) and darks is (r,c)
    else:
        d = data
    #use slices to map super-pixels into last component -- (rows x columns) -> (rows/2 x cols/2 x 4)
    d_stack = np.concatenate([ d[s][...,None] for s in slices ], axis=-1)
    #matrix multiply last dimensions of gains and data
    if gains is not None:
        cal_stack = np.einsum('...ij,...j', gains, d_stack) 
    else:
        cal_stack = d_stack
    del d_stack
    
    ## Interpolation ##
    i_stack = interpolate(cal_stack,slices,imethod)
    del cal_stack
    
    ## Compute Stokes Parameters ##
    #split stack into components so that we can use them in numexpr
    i0 = i_stack[...,0]
    i90 = i_stack[...,1]
    i45 = i_stack[...,2]
    i135 = i_stack[...,3]
    del i_stack
    
    #pre-allocate an array to hold stokes vectors
    s = np.ndarray(i0.shape + (3,))
    s[...,0] = ne.evaluate('(i0+i90+i45+i135)/2')
    s[...,1] = ne.evaluate('i0-i90')
    s[...,2] = ne.evaluate('i45-i135')
    
    return s

RawData = namedtuple('RawData',['roll','pitch','heading','times','exps','data','cuts'])

def load_file(fname):
    with h5py.File(fname,'r') as f:
        roll = array(f['compass/roll'])
        pitch = array(f['compass/pitch'])
        heading = array(f['compass/heading'])
        times = array(f['timestamps'])
        exps = array(f['exposures'])
        data = squeeze(array(f['data']))
        cuts = [0]
        return RawData(roll,pitch,heading,times,exps,data,cuts)

def load_file_shapes(fname):
    with h5py.File(fname,'r') as f:
        roll = f['compass/roll'].shape
        pitch = f['compass/pitch'].shape
        heading = f['compass/heading'].shape
        times = f['timestamps'].shape
        exps = f['exposures'].shape
        data = array(f['data'].shape)
        data = tuple(data[data>1]) #squeeze
        cuts = (1,)
        return RawData(roll,pitch,heading,times,exps,data,cuts)

def concatenate_shapes(ss):
    s1,s2 = ss
    if len(s1) == len(s2) and all([x==y for x,y in zip(s1[1:],s2[1:])]):
        return (s1[0]+s2[0],) + s1[1:]
    else:
        raise ValueError('Input array dimensions (except first) must match')

def load_dir(dname):
    files = glob(os.path.join(dname,'*.h5'))
    d = None
    i = 0
    cuts = []
    for f in files:
        print('{}... '.format(f),end='',flush=True)
        x = load_file(f)
        cuts.extend(i+c for c in x.cuts)
        i += x.times.shape[0]
        if d is None:
            d = x
        else:
            d = RawData._make(map(concatenate,zip(d,x)))
    return d._replace(cuts=cuts)

def process_data(rawdata,outfname,calfile=None,gps=None,caldata=None):
    
    if caldata is None:
        caldata = load_cal(calfile) #handles None case internally

    tmp_path = outfname+'.tmp'
    with h5py.File(tmp_path,'w') as outfile:
        outfile.create_dataset('roll',data=rawdata.roll)
        outfile.create_dataset('pitch',data=rawdata.pitch)
        outfile.create_dataset('heading',data=rawdata.heading)
        outfile.create_dataset('exposures',data=rawdata.exps)
        outfile.create_dataset('cuts',data=rawdata.cuts)
        if gps is not None:
            outfile.create_dataset('gps',data=gps)
        outfile.create_dataset('timestamps',data=rawdata.times)
        sout, dout, aout = None, None, None
        for i,frame in enumerate(rawdata.data):
            #process frame
            s = process(frame,caldata)
            if sout is None:
                sout = outfile.create_dataset('stokes',shape=(rawdata.data.shape[0],)+s.shape,dtype=float32,chunks=(1,)+s.shape,compression=5,shuffle=True)
            sout[i] = s
    if os.path.exists(outfname):
        os.rename(outfname,outfname+'.old')
        os.rename(tmp_path,outfname)
        os.remove(outfname+'.old')
    else:
        os.rename(tmp_path,outfname)

def process_file(fname,outfname=None,calfile=None,gps=None,caldata=None):
    try:
        rawdata = load_file(fname)
        
        if outfname is None:
            outfname = os.path.splitext(fname)[0]+'-stokes.h5'

        process_data(rawdata,outfname,calfile,gps,caldata)

    except KeyboardInterrupt as ki:
        raise ki
    except Exception as e:
        print('Failed:',e)

def process_dir(indir,outdir=None,calfile=None,gps=None,caldata=None,force=False):
    try:
        #load calibration data (handles None case properly)
        if caldata is None:
            caldata = load_cal(calfile)
        #all h5 files:
        files = glob(os.path.join(indir,'*.h5'))
        #check if output directory exists
        if outdir is None:
            outdir = indir
        elif not (os.path.exists(outdir) and os.path.isdir(outdir)):
            os.mkdir(outdir)
        #process each file in turn, as necessary
        for f in files:
            #input filename = 'xxx.h5' --> output to 'xxx-stokes.h5'
            outfname = os.path.join(outdir,os.path.basename(f)[:-3]+'-stokes.h5')
            if not force and os.path.exists(outfname) and file_newer(outfname,f):
                print('Already processed {}'.format(f))
                continue

            print('Processing {}...'.format(f),flush=True)
            process_file(f,outfname,gps=gps,caldata=caldata)            
            
    except KeyboardInterrupt as ki:
        raise ki
    except Exception as e:
        print('Failed:',e)
