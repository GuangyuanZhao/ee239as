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
import matplotlib as mpl
from bisect import bisect_left, bisect_right
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from pylab import *


#~ colormaps
def apply_cmap(X,cmap=None,norm=None,vmin=None,vmax=None,alpha=None,bytes=False):
    if norm is None:
        norm = Normalize(vmin,vmax)
    return cmap(norm(X),alpha,bytes)

def make_colormap(name,colors,N=None,gamma=1.0,typ='smooth'):
    """make a matplotlib colormap typ is either smooth or segmented"""
    if N is None:
        N = 256
    if typ == 'smooth':
        try:
            colors['red']
            cmap = LinearSegmentedColormap(name, colors,N,gamma)
        except:
            cmap = LinearSegmentedColormap.from_list(name,colors,N,gamma)
    elif typ == 'segmented':
        cmap = ListedColormap(colors,name,N)
    cm.register_cmap(name,cmap)
    return cmap

#red, orange, green, blue, purple; with subtly increasing lightness (chroma=1.08)
rbow1 = make_colormap('rbow1',['#BB654C','#B27A37','#999136','#75A552','#47B481','#1ABEB7','#59C1E4','#AFBBF8','#FDAEEE'],gamma=1.6)
#as above, but constant lightness, chroma=1.1
rbow2 = make_colormap('rbow2',['#EC7F60','#D19040','#A7A13B','#75AC57','#3CB183','#09B0AF','#55A8CC','#9D99D0','#D386B8'],gamma=1.6)
#red, purple, blue, green; with increasing lightness, chroma=1.08
rdgn = make_colormap('rdgn',['#A15246','#B55C6B','#B96F95','#A989BD','#85A6D8','#50C2E1','#1FDAD5','#4FEDB8','#95FC93'],gamma=1.6)
gnrd = make_colormap('gnrd',['#95FC93','#4FEDB8','#1FDAD5','#50C2E1','#85A6D8','#A989BD','#B96F95','#B55C6B','#A15246'],gamma=1.6)
#hue wheel, constant lightness, chroma=1.08
ihsv = make_colormap('ihsv',['#FD8A69','#D7A043','#9CB34C','#58BC7D','#19BCB6','#68B1DA','#BE9BD5','#F486AA','#FD8A69'],gamma=1.6)
#;

def _smooth_cmap(a,b,c,f=None):
    #smoothly vary colors from a,b,c using f
    if f is None:
        f = lambda x: sqrt(1+3*x**2)-1
    ra,ga,ba = cm.colors.colorConverter.to_rgb(a)
    rb,gb,bb = cm.colors.colorConverter.to_rgb(b)
    rc,gc,bc = cm.colors.colorConverter.to_rgb(c)
    def cf(a,b,c):
        return lambda x: where(x<0.5, b+(a-b)*f(1-2*x), b+(c-b)*f(2*x-1))
    return {'red': cf(ra,rb,rc),'green':cf(ga,gb,gc),'blue':cf(ba,bb,bc)}
f_cusp = lambda x: sqrt(((x+1)**2-1)/3)
RdKBu = make_colormap('RdKBu',_smooth_cmap('#ff3030','0','#3030ff',f_cusp),512)
RdKBuL = make_colormap('RdKBuL',_smooth_cmap('#ff8888','0.5','#8888ff',f_cusp),512)
BuKRd = make_colormap('BuKRd',_smooth_cmap('#3030ff','0','#ff3030',f_cusp),512)
BuKRdL = make_colormap('BuKRdL',_smooth_cmap('#8888ff','0.5','#ff8888',f_cusp),512)
OrKCy = make_colormap('OrKCy',_smooth_cmap('#ffa632','0','#32d3ff',f_cusp),512)
CyKOr = make_colormap('CyKOr',_smooth_cmap('#32d3ff','0','#ffa632',f_cusp),512)

#RdPBu = make_colormap('RdPBu',_smooth_cmap('#fb7565','#d482c3','#2aa9d7',f_cusp),512)
#RdPBu = make_colormap('RdPBu',_smooth_cmap('#f66157','#ce73c3','#149ee1',f_cusp),512)
RdPBu = make_colormap('RdPBu',_smooth_cmap('#f14245','#ce54b4','#1f86e0',f_cusp),512)

#adapted from http://nbviewer.jupyter.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
def colorline(x, y, z=None, cmap=None, vmin=None, vmax=None, scalex=True,scaley=True,**args):
    '''
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    '''
    cmap = cm.get_cmap(cmap)
    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x)-1)
        if vmin is None:
            vmin = 0.0
        if vmax is None:
            vmax = 1.0

    if vmin is None:
        vmin = np.nanmin(z)
    if vmax is None:
        vmax = np.nanmax(z)
    norm = Normalize(vmin,vmax)
    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])
        
    z = np.asarray(z)
    
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, array=z, cmap=cmap, norm=norm, **args)
    ax = gca()
    ax.add_collection(lc)
    ax.autoscale_view(scalex=scalex,scaley=scaley)
    return lc

def symlog(v, linthresh=1, linscale=1, base=10):
    """apply symlog transform to v
    v: array_like
    linthresh: for |v| > linthresh, use log scaling
    linscale: scale factor for |v| <= linthresh
    base: base of the logarithm

    adapted from matplotlib/scale.py
    """
    v = asarray(v)
    ls = linscale
    lt = linthresh
    lb = log(base)

    s = np.sign(v)
    m = ma.masked_inside(v, -lt, lt, copy=False)
    lg = s * lt * (ls + ma.log(np.abs(m) / lt) / lb)
    if m.mask.any():
        return ma.where(m.mask, v * ls, lg)
    else:
        return lg


def isymlog(v, linthresh=1, linscale=1, base=10):
    """apply inverse symlog transform to v"""
    v = asarray(v)
    ls = linscale
    lt = linthresh
    ilt = lt*ls #symlog(lt, lt, linscale, base)
    lb = log(base)

    s = np.sign(v)
    m = ma.masked_inside(v, -ilt, ilt, copy=False)
    ep = s * lt * (ma.power(base, (s * (m / lt)) - ls))
    if m.mask.any():
        return ma.where(m.mask, v / ls, ep)
    else:
        return ep

def blend(c1,c2,a,g=1.0,n=3):
    """blend(c1,c2,a,g=1.0,n=3)
    return the color of c1 and c2 blended together with weight a and gamma g
    if g == 1.0, the result is (1-a)*c1 + a*c2
    otherwise, the result is ((1-a)*(c1**g) + a*(c2**g))**(1/g)
    Only the first n components of the last dimension of c1 and c2 are treated
    components beyond n are copied from c1
    """
    c1,c2 = broadcast_arrays(c1,c2,subok=True)
    y = empty_like(c1)
    if c1.shape[-1] > n:
        y.T[n:] = c1.T[n:]
    if g == 1.0:
        y.T[:n] = (1-a)*c1.T[:n] + a*c2.T[:n]
    else:
        y.T[:n] = power((1-a)*power(c1.T[:n],g) + a*power(c2.T[:n],g), 1.0/g)
    return y

def blend_alpha(c1,c2,axis=-1,alpha_channel=3):
    """c1 over c2"""
    c1,c2 = broadcast_arrays(c1,c2,subok=True)
    c = empty_like(c1)
    c1_view = moveaxis(c1,axis,0)
    c2_view = moveaxis(c2,axis,0)
    c_view = moveaxis(c,axis,0)
    a1,a2 = c1_view[alpha_channel],c2_view[alpha_channel]
    a2 = a2*(1-a1)
    a = a1+a2
    c_view[alpha_channel] = a
    if isscalar(a):
        if a != 0:
            for i in range(c1_view.shape[0]):
                if i == alpha_channel: continue
                c_view[i] = (c1_view[i]*a1 + c2_view[i]*a2) / a
        else:
            c_view[:] = 0
    else:
        for i in range(c1_view.shape[0]):
            if i == alpha_channel: continue
            c_view[i] = (c1_view[i]*a1 + c2_view[i]*a2) / a
            c_view[i][a==0] = 0
    return c

def desat(c,a=0.5,n=3):
    """desat(c,a=0.5,n=3)
    desaturate colors by moving the first n color components by a towards their mean
    color components are assumed to be interleaved along the last axis of c
    """
    c = asarray(c)
    x = mean(c.T[:n],axis=0,keepdims=True)
    return blend(c,x,a,n=n)

def darken(c,a=0.5,n=3):
    """darken(c,a=0.5,n=3)
    darken a color by moving the first n color components towards 0
    color components are assumed to be interleaved along the last axis of c
    """
    return blend(c,0,a,n=n)

def lighten(c,a=0.5,w=1,n=3):
    """lighten(c,a=0.5,n=3)
    lighten a color by moving the first n color components towards w
    color components are assumed to be interleaved along the last axis of c
    """
    return blend(c,w,a,n=n)

def set_alpha(c,v,n=4):
    c = asarray(c)
    if c.shape[-1] < n:
        y = zeros(c.shape[:-1]+(n,))
        y.T[:c.shape[-1]] = c.T
    else:
        y = copy(c)
    y.T[n-1] = v
    return y

def ellipse(xx,yy,n=50,ax=None,xy0=(0,0),tup=True,**plot_args):
    """xx is major axis vector, yy is minor axis vector, n is points"""
    xx,yy = array(xx), array(yy)
    xf = concatenate((xx[None],yy[None])).T
    t = linspace(0,2*pi,n)
    ell = dot(xf,array((cos(t),sin(t))))
    if ax is not None:
        ax.plot(ell[0]+xy0[0],ell[1]+xy0[1],**plot_args)
    if not tup:
        return ell
    else:
        return (ell[0],ell[1])

def make_colorbar(cmap,vmin=None,vmax=None,norm=None,**kw):
    """make_colorbar(cmap,vmin=None,vmax=None,norm=None,**kw)
    creates a colorbar from a colormap and limits
    norm is None or a Normalize (or sub-class) object
      see matplotlib.colors.Normalize, PowerNorm, LogNorm, etc.
    if vmin or vmax is None, they will be taken from norm.
    if norm is also None, they will be 0 and 1, respectively
    additional arguments are forwarded to pyplot.colorbar()
    """
    if norm is None:
        if vmin is None:
            vmin = 0.0
        if vmax is None:
            vmax = 1.0
        norm = Normalize(vmin,vmax)
    if vmin is None:
        vmin = norm.vmin
    if vmax is None:
        vmax = norm.vmax
    sm = ScalarMappable(norm,cmap=cmap)
    sm.set_array([vmin,vmax])
    return colorbar(sm,**kw)

def make_colordisk(cmap,size,periods=1):
    pass
    """ mCDPeriods = periods;
        mCDisk = QImage(size, size, QImage::Format_ARGB32);
        //fill disk with colors
        for (int r = 0; r < size; ++r) {
            for (int c = 0; c < size; ++c) {
                // position 0 is at c == (size-1)/2
                double x = c - double(size - 1) / 2.0;
                double y = double(size - 1) / 2.0 - r;
                double d = hypot(x, y); //distance from center
                QRgb color;
                if (d < (size - 1)/2) {
                    double a = atan2(y, x) / 6.283185307179586; //-.5 to .5

                    if (a < 0) a += 1; //0 to 1
                    double ap; //will hold the period that a falls in
                    a = modf(a*periods, &ap);
                    int i = int(a*nextafter(mTableSize, 0)); //color index
                     color = mColorBuffer[i];
                }
                else {
                    color = qRgba(0, 0, 0, 0); //transparent black
                }
                *((QRgb*)mCDisk.bits() + r*size + c) = color;
                
            }
        }
        //make outline
        QPainter p(&mCDisk);
        p.setRenderHint(QPainter::Antialiasing);
        p.setBrush(Qt::NoBrush);
        qreal psz = size / 50.0;
        p.setPen(QPen(QBrush(Qt::black), psz));
        double r = size / 2.0; 
        QPointF center(r,r);
        p.drawEllipse(center, r-psz/2,r-psz/2);
        p.setBrush(QBrush(Qt::black));
        p.drawEllipse(center, r / 5, r / 5);
        p.end();"""

def align_labels(axes_list,axis='y',align=None):
    """align_labels(axes_list,axis='y',align=None)
    align either the x or y labels of all of the axes in axes_list
    axis must be either 'x' or 'y'
    align must be one of 'left','right','top','bottom','l','r','t','b'
    align defaults to 'left' if axis is 'y' and 'bottom' if axis is 'x'
    """
    if axis not in ('x','y'):
        raise ValueError("axis must be 'x' or 'y'")
    if align is None:
        align = 'l' if axis == 'y' else 'b'
    elif align not in ('l','r','t','b','left','right','top','bottom'):
        raise ValueError("invalid value for align argument")
    yx,xy = [],[]
    for ax in axes_list:
        yx.append(ax.yaxis.label.get_position()[0])
        xy.append(ax.xaxis.label.get_position()[1])
    
    if axis == 'x':
        if align in ('t','top'):
            lim = max(xy)
        elif align in ('b','bottom'):
            lim = min(xy)
    else:
        if align in ('l','left'):
            lim = min(yx)
        elif align in ('r','right'):
            lim = max(yx)

    if align in ('t','b','top','bottom'):
        for ax in axes_list:
            t = ax.xaxis.label.get_transform()
            x,y = ax.xaxis.label.get_position()
            ax.xaxis.set_label_coords(x,lim,t)
    else:
        for ax in axes_list:
            t = ax.yaxis.label.get_transform()
            x,y = ax.yaxis.label.get_position()
            ax.yaxis.set_label_coords(lim,y,t)

def legend_handle(label,**kw):
    return mpl.lines.Line2D([],[],label=label,**kw)

def text_bbox(text):
    r = text.figure.canvas.get_renderer()
    t = text.axes.transData.inverted()
    return t.transform_bbox(text.get_window_extent(r))


def tweak_x(x,y,w,h,eps=None):
    """return x tweaked so that points don't overlap
    x and y are iterables of the same length
    w,h are scalars, width and height, respectively, of points in data space
    eps is a tolerance for floating point comparisons. if None, defaults to machine epsilon
    returns a list of new x coordinates
    """
    if eps is None:
        wx2 = nextafter(w*2,0)
    else:
        wx2 = w*2-eps
    new_x = None
    x_y, y_y = None,None #points sorted by y values
    for px,py in zip(x,y):
        if new_x is None: #1st iteration
            new_x, x_y, y_y = [px],[px],[py]
            continue
        #find i,j such that for all y in y_y[i:j], py-h <= y <= py+h
        if eps is None:
            pymh,pyph = nextafter(py-h,inf), nextafter(py+h,-inf)
        else:
            pymh,pyph = py-h+eps, py+h-eps
        i = bisect_left(y_y,pymh)
        j = bisect_right(y_y,pyph,i)
        #special cases:
        if i == len(y_y) or j == 0: #all are < py-h or all are > py+h
            new_x.append(px); x_y.insert(i,px); y_y.insert(i,py)
            continue
        #sort by x, decreasing so that we prefer positive tweaks
        x_x = [inf] + sorted(x_y[i:j],reverse=True)+[-inf]
        #iterate thru candidate positions
        #and choose the closest
        nx = inf
        cps = []
        for k in range(len(x_x)-1):
            #is the gap between these points big enough?
            if x_x[k]-x_x[k+1] >= wx2:
                #is px in the gap?
                if x_x[k]-w >= px and px >= x_x[k+1]+w:
                    nx = px
                    break
                #left edge of the gap
                cp = x_x[k+1]+w
                if abs(cp-px) < abs(nx-px):
                    nx = cp
                #right edge of the gap
                cp = x_x[k]-w
                if abs(cp-px) < abs(nx-px):
                    nx = cp
        #insert
        k = bisect_left(y_y,py,i,j)
        new_x.append(nx)
        x_y.insert(k,nx)
        y_y.insert(k,py)
    return new_x

def tickparams(*args,**kw):
    ax = kw.pop('axes',None)
    if ax is None:
        ax = gca()
    return ax.tick_params(*args,**kw)

def _mticks(axis,locs=None,labels=None,**kw):
    tck = matplotlib.ticker
    class StrFmt(tck.Formatter):
        """As matplotlib.ticker.StrMethodFormatter, but does not require the format field to be labeled x"""
        def __init__(self,fmt):
            self.fmt = fmt
        def __call__(self,x,pos=None):
            return self.fmt.format(x)
    class FuncFmt(tck.Formatter):
        """Calls given function to format tick labels"""
        def __init__(self,fmt):
            self.fmt = fmt
        def __call__(self,x,pos=None):
            return self.fmt(x)

    nticks = kw.pop('nticks',None)
    tick1On = kw.pop('tick1On',None)
    tick2On = kw.pop('tick2On',None)
    minor = kw.pop('minor',False)
    reset = kw.pop('reset',False)

    if reset:
        if minor:
            axis.set_minor_locator(tck.NullLocator())
            axis.set_minor_formatter(tck.NullFormatter())
            axis.set_tick_params('minor',True)
            axis.isDefault_minloc = True
            axis.isDefault_minfmt = True
        else:
            axis.set_major_locator(tck.AutoLocator())
            axis.set_major_formatter(tck.ScalarFormatter())
            axis.set_tick_params('major',True)
            axis.isDefault_majloc = True
            axis.isDefault_majfmt = True
        axis.reset_ticks()
    
    if nticks is not None:
        #evenly divide axis by number of ticks
        #this does not play well with nonlinear scales!
        a,b = axis.get_view_interval()
        locs = linspace(a,b,nticks)

    if locs is not None:
        if isinstance(locs, tck.Locator):
            if minor: axis.set_minor_locator(locs)
            else: axis.set_major_locator(locs)
        else:
            if minor: axis.set_minor_locator(tck.FixedLocator(locs))
            else: axis.set_major_locator(tck.FixedLocator(locs))
    if minor: locs = axis.minor.locator()
    else: locs = axis.major.locator()

    txts = None
    if labels is not None:
        if isinstance(labels, tck.Formatter):
            if minor: axis.set_minor_formatter(labels)
            else: axis.set_major_formatter(labels)
        elif hasattr(labels,'format') and callable(labels.format):
            #use object with format method on each label
            if minor: axis.set_minor_formatter(StrFmt(labels))
            else: axis.set_major_formatter(StrFmt(labels))
        elif callable(labels):
            #call function for each label
            if minor: axis.set_minor_formatter(FuncFmt(labels))
            else: axis.set_major_formatter(FuncFmt(labels))
        else:
            #this makes a FixedFormatter and sets text properties
            txts = axis.set_ticklabels(labels,minor=minor,**kw)
    if txts is None:
        if minor:
            fmt = axis.minor.formatter
            ticks = axis.get_minor_ticks(len(locs))
        else:
            fmt = axis.major.formatter
            ticks = axis.get_major_ticks(len(locs))
        fmt.set_locs(locs)

        txt1, txt2 = [],[]
        for i, tick in enumerate(ticks):
            if tick1On is not None: tick.tick1On = tick1On
            if tick2On is not None: tick.tick2On = tick2On
            if tick.label1On:
                tick.label1.set_text(fmt(locs[i],i))
                tick.label1.update(kw)
                txt1.append(tick.label1)
            if tick.label2On:
                tick.label2.set_text(fmt(locs[i],i))
                tick.label2.update(kw)
                txt2.append(tick.label2)
        txts = txt1+txt2

    return locs, txts

def mxticks(*args,**kw):
    """mxticks(*args,**kw)
    Get or set the current x-axis tick locations and labels

    locs, labels = mxticks()
    mxticks(arange(6))
    mxticks(arange(3),['Tom','Dick','Harry'])

    Labels can be specified with a format string:
      mxticks(linspace(0,360,5),'{:.2e}°')
    
    Labels can be specified with a callable:
      mxticks(linspace(0.0,1.0,5),lambda x: '{:.0f}%'.format(x*100))

    Format current ticks:
      mxticks(labels='{:.2e}°')

    Locations may be specified with a matplotlib.ticker.Locator object
    Labels may be specified with a matplotlib.ticker.Formatter object
    To operate on a specific set of axes:
      ax = subplot(211)
      mxticks(..., axes = ax, ...)
    To operate on the minor ticks:
      mxticks(..., minor = True, ...)
    To reset the ticks, labels, and text properties before applying new settings:
      mxticks(..., reset = True, ...)
    Other keyword args are matplotlib.text.Text properties
    """
    if len(args) > 2:
        raise TypeError('Illegal number of arguments to mxticks')
    ax = kw.pop('axes',None)
    if ax is None:
        ax = gca()
    locs,labels = _mticks(ax.xaxis,*args,**kw)
    draw_if_interactive()
    return locs, silent_list('Text xticklabel',labels)

def myticks(*args,**kw):
    """myticks(*args,**kw)
    Get or set the current y-axis tick locations and labels

    locs, labels = myticks()
    myticks(arange(6))
    myticks(arange(3),['Tom','Dick','Harry'])

    Labels can be specified with a format string:
      myticks(linspace(0,360,5),'{:.2e}°')

    Labels can be specified with a callable:
      myticks(linspace(0.0,1.0,5),lambda x: '{:.0f}%'.format(x*100))

    Format current ticks:
      myticks(labels='{:.2e}°')

    Locations may be specified with a matplotlib.ticker.Locator object
    Labels may be specified with a matplotlib.ticker.Formatter object
    To operate on a specific set of axes:
      ax = subplot(211)
      myticks(..., axes = ax, ...)
    To operate on the minor ticks:
      myticks(..., minor = True, ...)
    To reset the ticks, labels, and text properties before applying new settings:
      myticks(..., reset = True, ...)
    Other keyword args are matplotlib.text.Text properties
    """
    if len(args) > 2:
        raise TypeError('Illegal number of arguments to myticks')
    ax = kw.pop('axes',None)
    if ax is None:
        ax = gca()
    locs,labels = _mticks(ax.yaxis,*args,**kw)
    draw_if_interactive()
    return locs, silent_list('Text yticklabel',labels)

def xticks_fmt(fmt=None,xt=None,*args,**kw):
    """Set xticks to xt, with labels generated by
     calling fmt.format(x) on each x in xt.
     If fmt is None, remove the labels
     If xt is None, use the current tick locations
    """
    if xt is None:
        xt = xticks()[0]
    if fmt is None:
        xl = []
    else:
        xl = [fmt.format(x) for x in xt]
    return xticks(xt,xl,*args,**kw)

def yticks_fmt(fmt=None,yt=None,*args,**kw):
    """Set yticks to yt, with labels generated by
     calling fmt.format(y) on each y in yt.
     If fmt is None, remove the labels
     If yt is None, use the current tick locations
    """
    if yt is None:
        yt = yticks()[0]
    if fmt is None:
        yl = []
    else:
        yl = [fmt.format(y) for y in yt]
    return yticks(yt,yl,*args,**kw)
