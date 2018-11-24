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
from make_cal import *
infname = '../cal150082.h5'
outfname = '../cal150082_fe.h5'

#cropped image pixels:
r0, c0 = 112, 534
h, w = 856, 856

#load pixel parameters from original calibration file
with h5py.File(infname) as infile:
    old_A = array(infile['parameters/A'])
    old_d = array(infile['parameters/d'])
    old_pattern = array(infile['parameters/pattern'])
    old_r0, old_c0 = infile['parameters/rc0']

old_h, old_w = old_d.shape
roi_r0 = r0-old_r0
roi_r1 = roi_r0 + h
roi_c0 = c0-old_c0
roi_c1 = roi_c0 + w

if roi_r0 < 0 or roi_r1 > old_h or roi_c0 < 0 or roi_c1 > old_w:
    print('input calibration file does not cover new pixels!')

roi = (slice(roi_r0,roi_r1), slice(roi_c0, roi_c1))

A = old_A[roi]
d = old_d[roi]

#changes to the pattern
pattern = old_pattern
if roi_r0 % 2 == 1: #swap rows?
    pattern = pattern[[1,0]]
if roi_c0 % 2 == 1: #swap columns?
    pattern = pattern[:,[1,0]]

#compute calibration matrices
order = [0,90,45,135]
starts = {pattern[p]:p for p in [(0,0),(0,1),(1,0),(1,1)]}
#slices for vectorizing super pixels in A: first 2 dimensions are image dimensions, last dimension is individual pixel's A vector
A_slices = [(slice(starts[a][0],None,2), slice(starts[a][1],None,2), slice(None)) for a in order]

#what's the ideal A for those 4 pixels?
orad = deg2rad(array(order))
A_ideal = 0.5*array([ones_like(orad), cos(2*orad), sin(2*orad), zeros_like(orad)])

C = compute_calibration(A, A_slices, A_ideal)

#write to file
with h5py.File(outfname) as outfile:
    opts = {'chunks':True, 'compression':9, 'shuffle':True, 'fletcher32':True} #dataset options
    grp = outfile.create_group('parameters')
    h5pat = grp.create_dataset('pattern',data=pattern)
    h5rc0 = grp.create_dataset('rc0',data=(r0,c0))
    h5drk = grp.create_dataset('d',data=d,dtype=float32,**opts)
    grp.create_dataset('A',data=A,dtype=float32,**opts)

    grp = outfile.create_group('calibration')
    grp['darks'] = h5drk
    grp.create_dataset('gains',data=C,dtype=float32,**opts)
    grp['attributes/pattern'] = h5pat
    grp.create_dataset('attributes/y0',data=r0)
    grp.create_dataset('attributes/x0',data=c0)
