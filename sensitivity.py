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
import h5py
from glob import glob
from collections import namedtuple
from estimate_position import *
from polprocess import process_dir

from stats_utils import circle_mask
from plot_utils import *
import scipy.stats as stats
from simulate import oceanstokes
from sunpos import sunpos_az as sunpos


#sensitivity analysis:
# 1) process files: xxx.h5 -> stokes/xxx-stokes.h5
# 2) sample AoPs: stokes/xxx-stokes.h5 -> stokes/traces.h5
# 3) unwrap AoPs
# 4) clip good segments
# 5) for each good segment:
#   6) for each subsegment of duration T:
#     7) linear regression
#     8) noise = original - linear regression
#     9) M distance using noise & linear reg
#     10) M dist threshold for significance
#     11) how long does it take to hit threshold?
# 12) what's the average time to hit threshold?

def gps_dist(gps0,gps1):
    """gps_dist(gps0,gps1)
    Distance in m between gps coordinates
    gps0,gps1 : lat,lon pairs in degrees
    """
    return WGS84.Inverse(gps0[0],gps0[1],gps1[0],gps1[1])['s12']

def sample_around_circle(stokes, center, angles, radius, mask_radius,rotate_aop=True,ydown=True):
    mask_size = 2*mask_radius+2 #to deal w/ offsets in 0..1
    sy = 1
    if ydown: sy = -1
    xy = array(center)+radius*stack([cos(angles),sy*sin(angles)],-1)
    aops, dolps = [], []
    for p, a in zip(xy, angles):
        x,y = p
        r,c = int(y), int(x)
        xo, yo = x-c, y-r
        mask = circle_mask(mask_radius,mask_size,(yo,xo))[None,...,None]
        count = sum(mask)
        r0,c0 = r-mask_radius, c - mask_radius
        s = squeeze(apply_over_axes(sum,stokes[:,r0:r0+mask_size,c0:c0+mask_size]*mask,(1,2)))
        aop = Stokes.aop(s,-1)
        if rotate_aop:
            #the horizon is tangential to the circle
            aop = wrap(aop - a - pi/2, -pi/2, pi/2)
        dolp = Stokes.dolp(s,-1)
        aops.append(aop)
        dolps.append(dolp)
    aops = stack(aops,-1)
    dolps = stack(dolps,-1)
    return aops,dolps,xy

def color_masks(size,cmap,xy,angles,mask_radius):
    #assumes none of the color dots overlap!
    mask_size = 2*mask_radius+2 #as in sample_around_circle
    img = zeros(size+(4,)) #assume 2d size
    for t,_xy in zip(angles,xy):
        x,y = _xy
        color = cmap(t/(2*pi)) #rgba
        r,c = int(y), int(x)
        xo,yo = x-c, y-r
        mask = circle_mask(mask_radius,mask_size,(yo,xo))
        r0,c0 = r-mask_radius, c-mask_radius
        img_slice = img[r0:r0+mask_size,c0:c0+mask_size]
        #to get the edges of the dots right, we need to set color, then alpha
        #the image starts as all zero, so we can just add
        img_slice[...,:3] += (mask>0)[...,None]*color[:3] #use binary mask to set color
        img_slice[...,3] += mask #mask is alpha
        #img_slice[:] = blend_alpha(mask*color,img_slice)
    return img

def sample_traces(fname,angles,center,radius,mask_radius):
    #center = [214,214]
    #radius = 186
    #mask_radius = 8
    #load relevant data from file
    aops,dolps = [],[]
    with h5py.File(fname) as infile:
        cuts = list(infile['cuts'])
        exps = array(infile['exposures'])
        gps = array(infile['gps'])
        head = array(infile['heading'])
        pitch = array(infile['pitch'])
        roll = array(infile['roll'])
        timestamps = array(infile['timestamps'])
        #load in intervals that fit in memory
        for i,j in zip(cuts,cuts[1:]+[None]):
            if len(cuts) > 1: print('  Sampling {}:{}...'.format(i,j),flush=True)
            ss = array(infile['stokes'][i:j])
            a,d,xy = sample_around_circle(ss,center,angles,radius,mask_radius)
            aops.append(a)
            dolps.append(d)
    aops = concatenate(aops,axis=1)
    dolps = concatenate(dolps,axis=1)

    return gps,xy,cuts,exps,head,pitch,roll,timestamps,aops,dolps

def adjust_bad_frames(cuts, bad_frames_by_video):
    bad_frames = []
    for i,c,bfs in zip(range(len(cuts)),cuts,bad_frames_by_video):
        for a,b in bfs:
            #deal with negative and None indexes if we're not in the last cut
            #None means "through end of video" in this case
            if i < len(cuts) - 1:
                if a < 0: a = cuts[i+1] + a
                if b is None: b = cuts[i+1] - cuts[i]
                elif b < 0: b = cuts[i+1] + b - cuts[i]
            if a >= 0: a += c
            if b is not None and b >= 0: b += c
            if bad_frames and bad_frames[-1][1] == a:
                bad_frames[-1][1] = b
            else:
                bad_frames.append([a,b])
    return bad_frames

def find_good_frames(length,bad_frames):
    good_frames = [[0]]
    for a,b in bad_frames:
        if b is None: b = length
        if a == good_frames[-1][0]:
            #the last good-frame interval starts at the same frame as a bad-frame interval
            #so move the start of the good-frame interval to the end of the bad-frame interval
            good_frames[-1][0] = b
        else:
            #the next bad-frame interval starts after the good-frame interval
            #so end the good-frame interval at the start of the bad-frame interval
            good_frames[-1].append(a)
            #and start the next good interval at the end of the current bad interval
            good_frames.append([b])
    if good_frames[-1][0] < length: #close the last interval
        good_frames[-1].append(length)
    else: #exclude if interval starts at the end
        del good_frames[-1]
    return good_frames

def make_trace_file(indir, tfname, bad_angles, bad_frames_by_video, good_frames):
    mask_radius = 8
    #sample AoPs around circle
    th = arange(0,360,6)
    #for a,b in bad_angles: #do this in final stats rather than here
    #    th = compress(logical_or(th < a, th >= b), th)
    th = deg2rad(th)

    #load *-stokes.h5
    files = sorted(glob(os.path.join(indir,'*-stokes.h5')))
    sampled_data = []
    cuts = []
    n = 0
    for f in files:
        print('Loading {}...'.format(f),flush=True)
        sampled_data.append(sample_traces(f,th,[214,214],186,mask_radius))
        cuts.extend(n+c for c in sampled_data[-1][2])
        n += sampled_data[-1][3].shape[0]
    #we now have a list of tuples of arrays
    #these parameters are constant:
    gps,xy = sampled_data[0][:2]
    #cuts is sampled_data[x][2], which we dealt with above
    #concatenate all of the datasets
    sampled_data = tuple(map(concatenate,zip(*sampled_data)))
    exps,head,pitch,roll,timestamps,aops,dolps = sampled_data[3:]
    
    #if we don't already have a list of good frames, figure it out
    if good_frames is None:
        #if we don't have a list of bad frames, then no frames are bad
        if bad_frames_by_video is None:
            bad_frames_by_video = []
        bad_frames = adjust_bad_frames(cuts,bad_frames_by_video)
        good_frames = find_good_frames(exps.shape[0],bad_frames)

    #write to file        
    with h5py.File(tfname) as tfile:
        tfile.create_dataset('cuts',data=cuts)
        tfile.create_dataset('exposures',data=exps)
        tfile.create_dataset('gps',data=gps)
        tfile.create_dataset('heading',data=head)
        tfile.create_dataset('pitch',data=pitch)
        tfile.create_dataset('roll',data=roll)
        tfile.create_dataset('timestamps',data=timestamps)
        tfile.create_dataset('aops',data=aops)
        tfile.create_dataset('dolps',data=dolps)
        tfile.create_dataset('xy',data=xy)
        tfile.create_dataset('angles',data=th)
        tfile.create_dataset('mask_radius',data=mask_radius)
        tfile.create_dataset('good_frames',data=array(good_frames))
        tfile.create_dataset('bad_angles',data=array(bad_angles))

StokesData = namedtuple('StokesData',['roll','pitch','heading','times','exps','stokes','cuts'])

def load_stokes_frame(sdir,cuts,frame_num):
    #which video does the frame_num appear in?
    #use next() to get the index of the first cut after the frame_num
    # subtract 1 to get the segment that contains the frame number
    vidx = next((i for i,v in enumerate(cuts) if v > frame_num), len(cuts))-1
    #subtract the cut to get the index of the frame inside the file
    fidx = frame_num - cuts[vidx]
    #find the associated *-stokes.h5 file
    fname = os.path.join(sdir,'{:03d}-stokes.h5'.format(vidx))
    #load the frame info from that video
    with h5py.File(fname,'r') as f:
        exp = f['exposures'][fidx]
        head = f['heading'][fidx]
        pitch = f['pitch'][fidx]
        roll = f['roll'][fidx]
        stokes = array(f['stokes'][fidx])
        time = f['timestamps'][fidx]
    return StokesData(roll,pitch,head,time,exp,stokes,[0])

TraceData = namedtuple('TraceData',['cuts','exposures','gps','heading','pitch','roll','timestamps','aops','dolps','xy','angles','mask_radius','good_frames','bad_angles'])

def load_traces(tfname):
    with h5py.File(tfname,'r') as tfile:
        cuts = list(tfile['cuts'])
        exposures = array(tfile['exposures'])
        gps = array(tfile['gps'])
        heading = array(tfile['heading'])
        pitch = array(tfile['pitch'])
        roll = array(tfile['roll'])
        timestamps = array(tfile['timestamps'])
        aops = array(tfile['aops'])
        dolps = array(tfile['dolps'])
        xy = array(tfile['xy'])
        angles = array(tfile['angles'])
        mask_radius = array(tfile['mask_radius'])[()] #this is a silly way to read a scalar
        good_frames = array(tfile['good_frames'])
        bad_angles = array(tfile['bad_angles'])
    return TraceData(cuts,exposures,gps,heading,pitch,roll,timestamps,aops,dolps,xy,angles,mask_radius,good_frames,bad_angles)

if __name__ == '__main__':

    #video_dir = '../../Data/2016-12-20 Argentina/2016-12-25-14.20.27-UTC+0000'
    #calfile = '../cal150082_fe.h5'
    #gps = (-31.150901,-64.342749)
    #elev = 805 #meters above WGS84
    #bad_frames_by_video = None
    #good_frames = [[6010,9600]] #videos 5-7
    #bad_angles = [[260,266]] #arm that holds sun shade

    #video_dir = '../../Data/2016-12-20 Argentina/2017-01-03-17.24.50-UTC+0000'
    #calfile = '../cal150082_fe.h5'
    #gps = (-31.150901,-64.342749)
    #elev = 805 #meters above WGS84
    #bad_frames_by_video = None
    #good_frames = [[25300,33300]] #videos 20-?
    #bad_angles = [[260,266]] #arm that holds sun shade

    # windy:
    #video_dir = 'O:/mexico_pol/2016-12-28-21.25.06-UTC+0000'
    #calfile = '../cal150082_fe.h5'
    #gps =  (-31.150901,-64.342749)
    #elev = 805
    #bad_frames_by_video = None
    #good_frames = [[900,13222]] #skip beginning
    #bad_angles = [[55,90],[260,285]] #over exposed areas at top & bottom

    
    video_dir = 'O:/mexico_pol/2017-03-22-16.05.01-UTC+0000'
    calfile = '../cal150082_fe.h5'
    gps =  (20.472043, -87.257618)
    elev = 0
    bad_frames_by_video = None
    good_frames = [[0,10700]] #all
    bad_angles = [[260,266]] #arm

    #video_dir = 'O:/mexico_pol/2017-01-03-17.24.50-UTC+0000'
    #calfile = '../cal150082_fe.h5'
    #gps =  (-31.150901,-64.342749)
    #elev = 805
    #bad_frames_by_video = None
    #good_frames = [[24040,33656]] #section that meets assumptions
    #bad_angles = [[260,266]] #arm for sunshade


    #sample AoPs and make trace file
    stokes_dir = os.path.join(video_dir,'stokes')
    tfname = os.path.join(stokes_dir,'traces.h5')
    if not os.path.exists(tfname):
        # process files: xxx.h5 -> stokes/xxx-stokes.h5
        process_dir(video_dir,stokes_dir,calfile,gps=gps)
        make_trace_file(stokes_dir,tfname,bad_angles,bad_frames_by_video,good_frames)

    #load trace file and unpack data
    td = load_traces(tfname)
    
    for good_frames in td.good_frames:
        #select data from the current interval
        time_slice = slice(*good_frames)
        aops = td.aops[time_slice]
        dolps = td.dolps[time_slice]
        t = td.timestamps[time_slice]/1e3 #posix timestamp, seconds
        te = t - t[0] #ellapsed time, seconds
        tc = td.timestamps[td.cuts]/1e3 - t[0] #ellapsed time of video cuts
        #mask bad angles:
        angle_mask = ones_like(td.angles,bool)
        for a,b in deg2rad(td.bad_angles):
            angle_mask *= logical_or(td.angles < a, td.angles >= b)
        #keep only traces with DoLP > 5% and good angles
        #dolp_mask = amin(dolps,axis=0) > 0.05
        #angle_mask *= dolp_mask

        aops = aops[:,angle_mask]
        dolps = dolps[:,angle_mask]
        angles = td.angles[angle_mask]
        xy = td.xy[angle_mask]

        #unwrap angles so that there are no discontinuities
        aops_uw = angle_unwrap(aops,axis=0,period=pi)

        #do linear regression on unwrapped angles
        te_1 = stack((te,ones_like(te)),-1) #shape (T, 2)
        a_x0 = lstsq(te_1, aops_uw)[0] #shape (2, k)
        #separate noise
        aops_lin = dot(te_1,a_x0)
        aops_noise = aops_uw - aops_lin
        #noise variance
        aops_nvar = angle_var(aops_noise,axis=0,period=pi)
        #inverse covariance matrix
        icov = inv(diag(aops_nvar))
        #icov = inv(cov(aops_noise.T)) #how does it go with this?

        #mahalanobis distance = sqrt((x0-x1).T * icov * (x0-x1))
        #  we want the mdist of each point from the initial point
        #  because the points have a linear relationship, we can simplify
        #  x(t) = dot(t,a) + x0
        #  mdist(x0, x(t)) = t * sqrt(a.T * icov * a)
        md = sqrt(einsum('i,ij,j', a_x0[0], icov, a_x0[0]))

        #how long will it take the distance change?
        #Null hypothesis: x(t) is drawn from the distribution with mean x0
        #Test statistic: mdist
        #significance: alpha = 5%, 1%, 0.1%
        #Find x' such that Prob(x(t) > x') = alpha, if x(t) > x' then reject null hypothesis
        # Prob(x(t) > x') is the survival function, use inverse to find x'
        #mdist^2 is chi^2, so use sqrt of chi2 inverse survival function
        alpha = array([5,1,0.1])/100
        md_thresholds = sqrt(stats.chi2.isf(alpha, aops.shape[1]))
        #when does t*md == md_thresholds? t = md_thresholds/md
        sig_times = md_thresholds/md

        #what gps variation do we get from the sun movement over sig_time?
        #look at sun positions at mean time of experiment +/- sig_times
        t_mean= mean(t)
        t_sp = t_mean+concatenate(([0],sig_times))
        sps = sunpos(t_sp,td.gps[0],td.gps[1],elev,radians=True)
        
        #how much did the sun move over that time?
        arcdists = array([arcdist(sps[0],sp) for sp in sps[1:]])

        #how much did the earth rotate over that time, in longitude?
        dlong = concatenate(([0], sig_times*360/86400))
        #find the gps positions that correspond to those sun positions at mean(t)
        sp_err = lambda gps,sp: arcdist(sp,sunpos(t_mean,gps[0],gps[1],elev,radians=True),True)
        sig_gps = []
        dist_m = []
        for sp in sps:
            sig_gps.append(minimize(sp_err,td.gps,args=(sp,)).x)
            dist_m.append(gps_dist(td.gps,sig_gps[-1]))
        dist_m = array(dist_m)
        dist_mi = units.convert(dist_m,'m','mi')

        fig1,ax = subplots()
        for a, th in zip(aops_uw.T,angles):
            c = cm.hsv(th/(2*pi))
            ax.plot(te,rad2deg(a),color=c)
        ymin,ymax=ax.get_ylim()
        #99% line:
        ax.vlines(sig_times[1],ymin,ymax,linestyles='dotted')
        #lines between videos
        ax.vlines(tc,ymin,ymax,linestyles='dashed')
        for i,x in enumerate(tc):
            ax.text(x+20,0.9*ymax,'{:03d}'.format(i))
    
        ax.set_ylim(ymin,ymax)
        ax.set_xlim(0,te[-1])
        ax.set_xlabel('Elapsed Time (s)')
        ax.set_ylabel('Polarization Angle')
        myticks(nticks=5,labels='{:.0f}°',axes=ax)

        #image:
        try:
            frame = load_stokes_frame(stokes_dir,td.cuts,1203) #good_frames[0])
            #apply colormap to intensity image:
            int_img = apply_cmap(frame.stokes[...,0], cm.gray, vmin=0, vmax=8096)
            #put colored dots on top:
            int_img = blend_alpha(color_masks(int_img.shape[:2],cm.hsv,xy,angles,td.mask_radius),int_img)
            #mask out periphery
            int_img[...,-1] *= circle_mask(int_img.shape[0]/2,int_img.shape[:2],center=True)
            fig2,ax = subplots()
            ax.imshow(int_img)
        except Exception as e:
            pass
        



def apply_filter_angle(noisy_signal,filter,period=2*pi):
    signal = empty_like(noisy_signal)
    coef = 2*pi/period
    for i in range(noisy_signal.shape[0]):
        signal[i] = angle(convolve(exp(1j*coef*noisy_signal[i]),filter,'same'))/coef
    noise = noisy_signal - signal
    return signal,noise

def apply_filter(noisy_signal,filter):
    """apply_filter(noisy_signal,filter)
    
    Parameters
    ----------
    noisy_signal : array_like, shape (N, T)
        
    filter : array_like, 1D
        Description
    
    Returns
    -------
    signal, noise : ndarray, shape (N, T)
    """
    signal = empty_like(noisy_signal)
    for i in range(noisy_signal.shape[0]):
        signal[i] = convolve(noisy_signal[i],filter,'same')
    noise = noisy_signal - signal
    return signal,noise

def apply_linreg_angle(noisy_signal,period=2*pi):
    """apply linear regression to angular data
    
    Parameters
    ----------
    noisy_signal : array_like, shape (N, T)
        regression is applied along T dimension
    period : scalar, optional (default 2*pi)
        period of angles. 2*pi for radians, 360 for degrees
    
    Returns
    -------
    signal, noise : array_like, shape (N,T)
    """
    #lstsq solves a x = b for x
    # a is input (M,N), b is output (M, K), x is coefficients (N,K)
    #t is input (time, 2)
    #ns_uw is output (time, channels)
    #we want (2, channels)
    ns_uw = angle_unwrap(noisy_signal,period=period)
    t = stack((arange(ns_uw.shape[1]),ones(ns_uw.shape[1])),-1)
    x = lstsq(t,ns_uw.T)[0]
    
    s = empty_like(noisy_signal)
    s.T[:] = dot(t,x)
    n = ns_uw - s
    return s,n

def synth_stokes(gps, elev, times, headings, pitch, fisheye=True):

    #get sun positions for times:
    # 5° above horizon after sunrise at t = 1482659100.0
    # 5 ° above horizone before sunset at t = 1482706500.0
    lat,lon = gps
    sp = sunpos(times,lat,lon,elev,radians=True)
    sun_az,sun_zn = sp[...,0], sp[...,1]
    stokes = oceanstokes(sun_az, sun_zn, headings[:,None], pitch)
    aops = Stokes.aop(stokes,-1)
    if fisheye:
        aops -= headings[:,None]
    return aops

def mean_mahalanobis_dist(signal, noise):
    '''signal and noise are shape (n x time)'''
    #variance along time
    noise_var = var(noise,axis=1)
    #covariance
    noise_cov = diag(noise_var)
    #invert to get Mahalanobis tranform matrix
    noise_icov = inv(noise_cov)
    #make pairwise squared M dist matrix
    #dist[r,c] = distance from point r to point r+c+1 in time
    #ie. r corresponds to time of the starting point, and c to the time difference between the starting and ending points
    #dist is upper left triangular, the lower right is filled with NaN
    dist = empty((signal.shape[1]-1,)*2)*nan
    for i in range(signal.shape[1]-1):
        #difference between point i and each subsequent point
        d = signal[:,i,None] - signal[:,i+1:]
        #compute mahalanobis distance = sqrt(transpose(d) * icov * d)
        m = sqrt(einsum('i...,ij,j...',d,noise_icov,d))
        dist[i,0:len(m)] = m
    #take square root
    #take mean over r dimension -- ie. average distance for each time difference between points
    d_mean = nanmean(dist,axis=0)
    #d_count = sum(dist*0+1,axis=0) #count non-NaN entries
    return d_mean, noise_var

def mahalanobis_dist_angle(signal, noise, window, period=2*pi):
    '''signal and noise are shape (time, n)
    returns mdist which is shape (time-window, window-1)
    and noise_var which is shape (time-window, n)
    '''
    T, n = signal.shape
    w = window
    mdist = empty(shape=(T-w, w-1))
    noise_var = empty(shape=(T-w, n))

    for i in range(T-w):
        #slice starting at current time-step
        sig = signal[i:i+w]
        nse = noise[i:i+w]
        #compute noise variance over window
        nse_var = angle_var(nse,axis=0,period=period)
        noise_var[i] = nse_var
        #form covariance matrix and invert
        nse_cov = diag(nse_var)
        nse_icov = inv(nse_cov)
        #compute M distance from i to each point in the window
        #d is difference, m is sqrt(transpose(d) * icov * d)
        d = angle_diff(sig[0], sig[1:], period)
        m = sqrt(einsum('...i,ij,...j',d,nse_icov,d))
        mdist[i] = m
    return mdist,noise_var

def mean_mahalanobis_dist_angle(signal, noise,period=2*pi):
    '''signal and noise are shape (n x time)'''
    #variance along time
    noise_var = angle_var(noise,axis=1,period=period)
    #covariance
    noise_cov = diag(noise_var)
    #invert to get Mahalanobis tranform matrix
    noise_icov = inv(noise_cov)
    #make pairwise squared M dist matrix
    #dist[r,c] = distance from point r to point r+c+1 in time
    #ie. r corresponds to time of the starting point, and c to the time difference between the starting and ending points
    #dist is upper left triangular, the lower right is filled with NaN
    dist = empty((signal.shape[1]-1,)*2)*nan
    for i in range(signal.shape[1]-1):
        #difference between point i and each subsequent point
        d = angle_diff(signal[:,i,None], signal[:,i+1:], period)
        #compute mahalanobis distance = sqrt(transpose(d) * icov * d)
        m = sqrt(einsum('i...,ij,j...',d,noise_icov,d))
        dist[i,0:len(m)] = m
    #take square root
    #take mean over r dimension -- ie. average distance for each time difference between points
    d_mean = nanmean(dist,axis=0)
    #d_count = sum(dist*0+1,axis=0) #count non-NaN entries
    return d_mean, noise_var

#francesc ferrer

#videos = [
#    ('../../Data/2016-10-02 Mermet/2016-10-02-19.08.47-UTC+0000/stokes.h5', #filename
#        [[[0,851]],[[868,1111]],[[10,241],[455,717],[900,1014],[1155,1202]],[[0,16],[175,246]],[[254,341],[1190,1202]],[[0,96],[455,486],[533,592],[750,946],[1040,1068],[1160,1202]]], #bad frames
#        [[60,120],[240,300]] #bad angles, blooming
#    ),
#]


#video_dir = '../../Data/2016-12-20 Argentina/2016-12-25-14.20.27-UTC+0000/stokes'
#videos = [os.path.join(video_dir,f) for f in os.listdir(video_dir) if f.endswith('stokes.h5')]

#~def make_trace_file(fname, tfname, bad_frames_by_video, bad_angles):
#~    #sample AoPs around circle, avoid top and bottom because of bloom
#~    th = arange(0,360,6)
#~    for a,b in bad_angles:
#~        th = compress(logical_or(th < a, th >= b), th)
#~    #th = deg2rad(concatenate((arange(0,60,6),arange(120,240,6),arange(300,360,6))))
#~    th = deg2rad(th)
#~
#~    print(fname,flush=True)
#~    aops,dolps = [],[]
#~    with h5py.File(fname) as infile:
#~        cuts = list(infile['cuts'])
#~        exps = array(infile['exposures'])
#~        gps0 = array(infile['gps'])
#~        head = array(infile['heading'])
#~        pitch = array(infile['pitch'])
#~        roll = array(infile['roll'])
#~        timestamps = array(infile['timestamps'])
#~        #load in intervals that fit in memory
#~        for i,j in zip(cuts,cuts[1:]+[None]):
#~            print('Processing {}:{}...'.format(i,j),flush=True)
#~            ss = array(infile['stokes'][i:j])
#~            a,d,xy = sample_around_circle(ss,[214,214],th,186,8)
#~            aops.append(a)
#~            dolps.append(d)
#~    aops = concatenate(aops,axis=1)
#~    dolps = concatenate(dolps,axis=1)
#~
#~    tso = timestamps/1000
#~    
#~    bad_frames = adjust_bad_frames(cuts,bad_frames_by_video)
#~    good_frames = find_good_frames(tso.shape[0],bad_frames)
#~        
#~    with h5py.File(tfname) as tfile:
#~        tfile.create_dataset('cuts',data=cuts)
#~        tfile.create_dataset('exposures',data=exps)
#~        tfile.create_dataset('gps',data=gps0)
#~        tfile.create_dataset('heading',data=head)
#~        tfile.create_dataset('pitch',data=pitch)
#~        tfile.create_dataset('roll',data=roll)
#~        tfile.create_dataset('timestamps',data=timestamps)
#~        tfile.create_dataset('aops',data=aops)
#~        tfile.create_dataset('dolps',data=dolps)
#~        tfile.create_dataset('xy',data=xy)
#~        tfile.create_dataset('angles',data=th)
#~        tfile.create_dataset('good',data=array(good_frames))

def concat_trace_files(fnames,bad_frames_by_video):
        tfile.create_dataset('cuts',data=cuts)
        tfile.create_dataset('exposures',data=exps)
        tfile.create_dataset('gps',data=gps0)
        tfile.create_dataset('heading',data=head)
        tfile.create_dataset('pitch',data=pitch)
        tfile.create_dataset('roll',data=roll)
        tfile.create_dataset('timestamps',data=timestamps)
        tfile.create_dataset('aops',data=aops)
        tfile.create_dataset('dolps',data=dolps)
        tfile.create_dataset('xy',data=stack((xx,yy)))
        tfile.create_dataset('angles',data=th)
        tfile.create_dataset('good',data=array(good_frames))





if False:
    #make ###-stokes.h5 and ###-traces.h5 files
    files = [os.path.join(path,f) for f in os.listdir(path) if f.endswith('.h5')]
    stokes_files = []
    for fn in files:
        print(fn,flush=True)
        outfile = os.path.join(path,'stokes',os.path.splitext(os.path.basename(fn))[0]+'-stokes.h5')
        stokes_files.append(outfile)
        if os.path.exists(outfile):
            print('Already done',flush=True)
            continue
        rphdtec = load_file(fn)
        process_data(rphdtec,outfile,None,gps,None)
    for fn in stokes_files:
        print(fn,flush=True)
        p,f = os.path.split(fn)
        outfile = os.path.join(p,f[:3]+'-traces.h5')
        if os.path.exists(outfile):
            print('Already done',flush=True)
            continue
        make_trace_file(os.path.join(path,fn),outfile,[],[[260,266]])

if False:
    #fname = '../../Data/2016-12-20 Argentina/2016-12-25-14.20.27-UTC+0000/stokes/stokes.h5'
    #bad_frames_by_video = [[(0,330)]] + [[]]*25 + [[(600,None)],[(0,None)],[(0,None)]]
    #bad_angles = [[260,266]] #arm that holds sun shade
    #tfname = os.path.join(os.path.dirname(fname),'traces.h5')
    #make_trace_file(fname,tfname,bad_frames_by_video,bad_angles)

    fname = '../../Data/2016-12-20 Argentina/2016-12-25-14.20.27-UTC+0000/stokes/traces.h5'
    vname = '2016-12-25-14.20.27-UTC+0000'
    with h5py.File(fname,'r') as infile:
        tso = array(infile['timestamps'])/1000.0 #divide by 1000 to get seconds
        aops = array(infile['aops'])
        dolps = array(infile['dolps'])
        xx,yy = array(infile['xy'])
        th = array(infile['angles'])
        good_frames = array(infile['good'])
        gps0 = array(infile['gps'])
        cuts = array(infile['cuts'])


    good_frames = [[6010,9600]] #videos 5-7

if False:
    #filter aop traces to get noise-free signals
    filt = sinc_filter(25,3)
    aops_f, aops_n = apply_filter_angle(aops, filt, pi)
    dolps_f, dolps_n = apply_filter(dolps, filt)

    #select the first good time-slice
    time_slice = slice(*good_frames[0])
    ts = tso[time_slice] #timestamps
    te = ts-ts[0] #elapsed time
    aops = aops[:,time_slice]
    aops_f = aops_f[:,time_slice]
    aops_n = aops_n[:,time_slice]
    dolps = dolps[:,time_slice]
    dolps_f = dolps_f[:,time_slice]
    dolps_n = dolps_n[:,time_slice]

if False:
    video_dir = '../../Data/2016-12-20 Argentina/2016-12-25-14.20.27-UTC+0000/stokes'
    videos = [os.path.join(video_dir,f) for f in os.listdir(video_dir) if f.endswith('stokes.h5')]

    bad_frames_by_video = []
    bad_angles = [[260,266]] #arm that holds sun shade

    times_all = []
    angles_all = []
    aops_f_all = []
    aops_n_all = []
    aops_v_all = []
    dolps_f_all = []
    dolps_n_all = []
    dist_all = []
    prob_all = []

    for fname in videos[8:]:
        vname = os.path.basename(fname)
        print(fname,flush=True)
        with h5py.File(fname) as infile:
            tso = array(infile['timestamps'])/1000
            ss = array(infile['stokes'])
            cuts = list(infile['cuts'])
            gps0 = array(infile['gps'])
        
        bad_frames = adjust_bad_frames(cuts,bad_frames_by_video)
        good_frames = find_good_frames(tso.shape[0],bad_frames)

        #sample AoPs around circle, avoid top and bottom because of bloom
        th = arange(0,360,6)
        for a,b in bad_angles:
            th = compress(logical_or(th < a, th >= b), th)
        #th = deg2rad(concatenate((arange(0,60,6),arange(120,240,6),arange(300,360,6))))
        th = deg2rad(th)
        aops,dolps,xy = sample_around_circle(ss,[214,214],th,186,8)
        
        #filter aop traces to get noise-free signals
        filt = sinc_filter(25,3)
        aops_f, aops_n = apply_filter_angle(aops, filt, pi)
        dolps_f, dolps_n = apply_filter(dolps, filt)

        #select the first good time-slice
        time_slice = slice(*good_frames[0])
        ts = tso[time_slice] #timestamps
        te = ts-ts[0] #elapsed time
        aops = aops[:,time_slice]
        aops_f = aops_f[:,time_slice]
        aops_n = aops_n[:,time_slice]
        dolps = dolps[:,time_slice]
        dolps_f = dolps_f[:,time_slice]
        dolps_n = dolps_n[:,time_slice]

        #select traces with dolp > 5% during the good time slice
        space_slice = find(all(dolps_f>0.05,axis=1))
        th = th[space_slice]
        aops = aops[space_slice]
        aops_f = aops_f[space_slice]
        aops_n = aops_n[space_slice]
        dolps = dolps[space_slice]
        dolps_f = dolps_f[space_slice]
        dolps_n = dolps_n[space_slice]

        dist, aops_v = mean_mahalanobis_dist_angle(aops_f,aops_n,pi)

        #probability that aops have changed:
        prob = stats.chi2.cdf(dist**2,aops.shape[0])
        #to find M distances that correspond to probabilities use
        #sqrt(stats.chi2.isf(1-prob,aops.shape[0]))

        try:
            i_sig = []
            i_sig.append(find(prob>=0.9)[0])
            i_sig.append(find(prob>=0.99)[0])
            i_sig.append(find(prob>=0.999)[0])
        except:
            pass
        t_sig = te[i_sig]
        p_sig = prob[i_sig]
        d_sig = dist[i_sig]

        times_all.append(ts)
        angles_all.append(th)
        aops_f_all.append(aops_f)
        aops_n_all.append(aops_n)
        dolps_f_all.append(dolps_f)
        dolps_n_all.append(dolps_n)
        aops_v_all.append(aops_v)
        prob_all.append(prob)
        dist_all.append(dist)

        fig,ax = subplots(2,1,num=vname)
        tlim = prob.shape[0]-1
        ax[0].plot(te[:tlim],1-prob[:tlim],'k')
        ax[0].spines['left'].set_color('k')
        ax[0].tick_params(axis='y',colors='k')
        ax[0].set_yscale('symlog',linthreshy=1e-4)
        ax[0].set_ylabel('Prob AoP Unchanged',color='k')
        ax[0].set_xlim((0,te[tlim-1]))

        for t,p in zip(t_sig,p_sig):
            ax[0].plot([te[0],t],[1-p,1-p],':k')
            ax[0].axvline(t,ls=':',c='k')

        for a,av,t in zip(aops_f,aops_v,th):
            c = cm.hsv(t/(2*pi))
            ax[1].plot(te,rad2deg(a),color=c)
            ax[1].fill_between(te,rad2deg(a-1*sqrt(av)),rad2deg(a+1*sqrt(av)),facecolor=desat(c),edgecolor='None',alpha=0.25)
        ax[1].set_ylim((-90,90))
        ax[1].set_xlim((0,te[-1]))
        ax[1].set_ylabel('AoP')
        ax[1].set_xlabel('Time (s)')
        fig.tight_layout()
        fig.savefig(vname+'-prob.png',dpi=300)
        close(fig)

        #fig=figure(vname+'-img')
        #ss_img = apply_cmap(ss[0,:,:,0],cm.cubehelix)
        #ss_img = blend_alpha(color_masks(ss_img.shape[:2],cm.hsv,xy,th,8,18),ss_img)
        #imshow(ss_img)
        #fig.tight_layout()
        #fig.savefig(vname+'-img.png',dpi=300)
        continue

        try:
            i_sig = []
            i_sig.append(find(prob>=0.9)[0])
            i_sig.append(find(prob>=0.99)[0])
            i_sig.append(find(prob>=0.999)[0])
        except:
            pass
        t_sig = te[i_sig]
        p_sig = prob[i_sig]
        d_sig = dist[i_sig]
        #limits for 40d chi2:
        # 99% : d = 7.9807
        # 99.9% : d = 8.5675
        fig=figure(fname+'-1')
        ax=gca()
        tlim=-1
        ax.plot(te[:tlim],1-prob[:tlim],'k')
        ax.spines['left'].set_color('k')
        ax.tick_params(axis='y',colors='k')
        ylabel('Probability AoP is Unchanged',color='k')
        xlabel('Time (s)')

        for t,p in zip(t_sig,p_sig):
            plot([te[0],t],[1-p,1-p],':k')
        
        yscale('symlog',linthreshy=1e-4)
        ylim(1e-4,1.1)
        ax2=twinx()
        ax2.spines['left'].set_color('r')
        ax2.spines['right'].set_color('b')
        ax2.plot(te[:tlim],dist[:tlim],'b')
        ax2.tick_params(axis='y',colors='b')
        ylabel('Mean AoP Mahalanobis Distance',color='b')
        for t,d in zip(t_sig,d_sig):
            plot([t,te[-1]],[d,d],':k')
            axvline(t,ls=':',c='k')
        
        xlim(0,te[tlim-1])
        title('AoP Sensitivity ({} traces)'.format(len(space_slice)))


        #gps0 = (37.284538,-88.858304) #Mermet Springs
        sp = sunpos(ts,gps0[0],gps0[1],105,radians=True)
        #arcdistance between sp[i] and sunposition at time 0 with given coordinates
        sp_err = lambda gps,i: arcdist(sp[i],sunpos(ts[0],gps[0],gps[1],105,radians=True),True)
        gps_sig = []
        dist_m = []
        for i in i_sig:
            gps_sig.append(minimize(sp_error,gps0,args=(i)).x)
            dist_m.append(gps_dist(gps0,gps_sig[-1]))
        dist_m = array(dist_m)
        dist_mi = units.convert(dist_m,'m','mi')

        figure(fname+'-2')
        for a,v,t in zip(aops_f,aops_v,th):
            c = cm.hsv(t/(2*pi))
            plot(te,rad2deg(a),color=c)
            fill_between(te,rad2deg(a-1*sqrt(v)),rad2deg(a+1*sqrt(v)),facecolor=desat(c),edgecolor='None',alpha=0.25)
            #fill_between(te,rad2deg(a-2*sqrt(v)),rad2deg(a+2*sqrt(v)),facecolor=desat(c),edgecolor='None',alpha=0.25)
            #fill_between(te,rad2deg(a-3*sqrt(v)),rad2deg(a+3*sqrt(v)),facecolor=desat(c),edgecolor='None',alpha=0.25)
        ylim(-25,55)
        axvline(23.2,ls=':',c='k')
        axvline(26.5,ls=':',c='k')
        axvline(28.9,ls=':',c='k')
        xlim(0,40)
        ylabel('Angle of Polarization')
        xlabel('Time (s)')
        myticks(labels='{:.0f}\xB0')

        
        figure(fname+'-3')
        ss_img = apply_cmap(ss[0,:,:,0],cm.cubehelix)
        ss_img = blend_alpha(color_masks(ss_img.shape[:2],cm.hsv,xy,th,8,18),ss_img)
        imshow(ss_img)
