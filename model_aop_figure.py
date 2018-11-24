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
from estimate_position import *
from matplotlib import colors

def stokes_to_rgb(s,ss=0.9,vs=0.2):
    a = wrap(Stokes.aop(s,-1),0,pi)
    p = Stokes.dolp(s,-1)

    sm,vm = 1-ss,1-vs
    rgb = colors.hsv_to_rgb(stack((a/pi,ss*p+sm,vs*p+vm),axis=-1))
    return rgb,a,p

def stokes_disk(n,ss=0.9,vs=0.2,r=(0.1,0.95,1),aa=4):
    r0,r1,r2=r
    x = linspace(-1,1,n)
    x,y = meshgrid(x,-x)
    r = hypot(x,y)
    a = arctan2(y,x)
    a = wrap(a,0,pi)/pi #two periods, scale 0 to 1
    p = clamp((r-r0)/(r1-r0)) #DoLP, 0 at r0, 1 at r1
    #colors
    sm,vm = 1-ss,1-vs
    #keep alpha channel separate
    al = ones_like(x)
    rgb = colors.hsv_to_rgb(stack((a,ss*p+sm,vs*p+vm),axis=-1))
    #mask inside r0 and outside r1 to black, but keep alpha
    m = circle_mask(r1*n/2,n,center=True) - circle_mask(r0*n/2,n,center=True)
    rgb = rgb*m[...,None]
    #mask outside r=1 to transparent
    al = al*circle_mask(r2*n/2,n,center=True)
    rgba = concatenate((rgb,al[...,None]),-1)
    return rgba,a,p

print('Figure 1',flush=True)
heading = deg2rad(linspace(-180,180,361))
pitch = deg2rad(linspace(-45,45,91))
H,P = meshgrid(heading,pitch)

print(' stokes at 80',flush=True)
stokes80 = oceanstokes(0,deg2rad(90-80),H,P)
print(' stokes at 45',flush=True)
stokes45 = oceanstokes(0,deg2rad(90-45),H,P)
print(' stokes at 10',flush=True)
stokes10 = oceanstokes(0,deg2rad(90-10),H,P)

print(' to RGB',flush=True)
rgb80,aop80,dolp80 = stokes_to_rgb(stokes80,1,0)
rgb45,aop45,dolp45 = stokes_to_rgb(stokes45,1,0)
rgb10,aop10,dolp10 = stokes_to_rgb(stokes10,1,0)

quiver_opts={'pivot':'mid','scale_units':'x','scale':0.1,'headlength':0,'headwidth':1,'width':0.003}

print(' figure',flush=True)
fig1=figure()
subplot(311)
imshow(rgb80,origin='lower',extent=(-180,180,-45,45))
#plot([-180,180],[0,0],'k:',scalex=False,scaley=False)
#axvline(90,ls=':',c='k')
#axvline(0,ls=':',c='k')
#axvline(-90,ls=':',c='k')
quiver(rad2deg(heading[::45]),rad2deg(pitch[5::20]),cos(aop80[5::20,::45]),sin(aop80[5::20,::45]),**quiver_opts)
mxticks(linspace(-180,180,9),labels='')
myticks(linspace(-40,40,5),labels='{:.0f}\xB0')
tickparams(length=0)

subplot(312)
ylabel('Detector Pitch')
imshow(rgb45,origin='lower',extent=(-180,180,-45,45))
#plot([-180,180],[0,0],'k:',scalex=False,scaley=False)
#axvline(90,ls=':',c='k')
#axvline(0,ls=':',c='k')
#axvline(-90,ls=':',c='k')
quiver(rad2deg(heading[::45]),rad2deg(pitch[5::20]),cos(aop45[5::20,::45]),sin(aop45[5::20,::45]),**quiver_opts)
mxticks(linspace(-180,180,9),labels='')
myticks(linspace(-40,40,5),labels='{:.0f}\xB0')
tickparams(length=0)

subplot(313)
imshow(rgb10,origin='lower',extent=(-180,180,-45,45))
#plot([-180,180],[0,0],'k:',scalex=False,scaley=False)
#axvline(90,ls=':',c='k')
#axvline(0,ls=':',c='k')
#axvline(-90,ls=':',c='k')
quiver(rad2deg(heading[::45]),rad2deg(pitch[5::20]),cos(aop10[5::20,::45]),sin(aop10[5::20,::45]),**quiver_opts)
mxticks(linspace(-180,180,9),labels='{:.0f}\xB0')
myticks(linspace(-40,40,5),labels='{:.0f}\xB0')
tickparams(length=0)

xlabel('Heading relative to Sun')

## color scale



print('Figure 2',flush=True)
##left vs right eye
eye_heading = 40
stokes_r = oceanstokes(0,deg2rad(90-linspace(10,80,8))[:,None],heading+deg2rad(eye_heading),0)
stokes_l = oceanstokes(0,deg2rad(90-linspace(10,80,8))[:,None],heading-deg2rad(eye_heading),0)
aop_r = Stokes.aop(stokes_r,-1)
aop_l = Stokes.aop(stokes_l,-1)
dolp_r = Stokes.dolp(stokes_r,-1)
dolp_l = Stokes.dolp(stokes_l,-1)

fig2=figure()

for al,ar in zip(aop_l,aop_r):
    plot(rad2deg(heading),rad2deg(al),'r')
    plot(rad2deg(heading),rad2deg(ar),'b')

xlim(-180,180)
mxticks(linspace(-180,180,9),labels='{:.0f}\xB0')
ylabel('Polarization Angle')
myticks(labels='{:.0f}\xB0')
axvline(ls=':',c='k')
xlabel('Animal heading relative to sun')
title('Angular difference between eyes: {}\xB0'.format(2*eye_heading))
text(130,50,'Left Eye',ha='center')
text(45,50,'Right Eye',ha='center')
xlabel('$\phi$, animal heading relative to sun')

text(-90,0,'Sun elevation\n$\\alpha$=80°',ha='center')
text(-90,-50,'$\\alpha$=10°',ha='center')

print('Figure 3',flush=True)
fig3=figure()
for al,ar in zip(aop_l,aop_r):
    plot(rad2deg(al),rad2deg(ar))
ax=gca()
ax.set_aspect('equal')
plot(rad2deg(aop_l[:,15::15]),rad2deg(aop_r[:,15::15]),':k');
title('Heading difference between eyes: {}\xB0'.format(2*eye_heading))
xlabel('Left eye e-vector angle')
ylabel('Right eye e-vector angle')
mxticks(labels='{:.0f}\xB0')
myticks(labels='{:.0f}\xB0')


print('Figure 4',flush=True)
fig4 = figure()
aop_h = oceanaop(deg2rad(60),deg2rad(45),heading,0)
aop_h_tan = aop_h - heading
ax = fig4.add_subplot(111,projection='polar')
ax.set_theta_offset(pi/2) #N up
ax.set_theta_direction(-1) #+ is right

quiver(heading[::15],ones_like(heading[::15]),cos(aop_h_tan[::15]),sin(aop_h_tan[::15]),**quiver_opts)
ax.set_rlim(0,1.1)
ax.set_rticks([])

print('Figure 5',flush=True)
### Strip of AoPs ###
#the x axis is angle, y is linear
#projecting onto a cylinder, radius=100, height = 200

#editable parameters:
radius, height = 1, 2
quiver_dh, quiver_dy = 15, 0.2 #quiver spacings in heading (deg), height
spu = 100 #model samples per linear unit
pps = 4 #pixels per sample in final render
sun_elevation = 10
bounds = (-180,180) #in degrees, 0 is to the sun

#computed parameters:
width = deg2rad(bounds[1]-bounds[0])*radius #arc distance = angle*radius
head = deg2rad(linspace(bounds[0],bounds[1],round(width*spu)))
y = linspace(0,height,round(height*spu)) - height/2
pitch = arctan2(y,radius)

H,P = meshgrid(head,pitch)
stokes = oceanstokes(0,deg2rad(90-sun_elevation),H,P)
rgb,aop,dolp = stokes_to_rgb(stokes)

fig5 = figure(figsize=rgb.T.shape[1:],dpi=pps)
ax = plt.Axes(fig5,[0.,0.,1.,1.])
ax.set_axis_off()
fig5.add_axes(ax)
ax.imshow(rgb,origin='lower',aspect='auto',extent=(head[0],head[-1],y[0],y[-1]))
fig5.savefig('AoP Cylinder, {} to {}, sun {}.png'.format(bounds[0],bounds[1],sun_elevation),dpi=pps) 
close(fig5)

#we want quiver lines every 15 degrees in heading
# and every 0.2 units of height
head_q = deg2rad(arange(bounds[0],bounds[1]+1,quiver_dh))
y_q = arange(0, height/2, quiver_dy)
y_q = concatenate((-y_q[:0:-1],y_q))
pitch_q = arctan2(y_q,radius)
H_q,P_q = meshgrid(head_q,pitch_q)
aop = Stokes.aop(oceanstokes(0,deg2rad(90-sun_elevation),H_q,P_q),-1)

#make a figure that's precisely the right number of pixels
# increase dpi slightly to get some interpolation
quiver_opts={'pivot':'mid','units':'width','scale_units':'width','scale':1.5*360/15,'headlength':0,'headaxislength':0,'headwidth':0.001,'width':0.001}
fig5 = figure(figsize=rgb.T.shape[1:],dpi=pps)
ax = plt.Axes(fig5,[0.,0.,1.,1.])
ax.set_axis_off()
fig5.add_axes(ax)
ax.imshow(rgb,origin='lower',aspect='auto',extent=(head[0],head[-1],y[0],y[-1]))
ax.quiver(head_q,y_q,cos(aop),sin(aop),**quiver_opts)
fig5.savefig('AoP Cylinder, {} to {}, sun {}, quivers.png'.format(bounds[0],bounds[1],sun_elevation),dpi=pps) 
close(fig5)

#quiver_opts={'pivot':'mid','units':'width','scale_units':'width','scale':1.5*360/15,'headlength':0,'headaxislength':0,'headwidth':0.001,'width':0.001}
#fig5 = figure(figsize=rgb.T.shape[1:],dpi=pps)
#ax = plt.Axes(fig5,[0.,0.,1.,1.])
#ax.set_axis_off()
#fig5.add_axes(ax)
#ax.imshow(ones_like(rgb),origin='lower',aspect='auto',extent=(head[0],head[-1],y[0],y[-1]))
#ax.quiver(head_q,y_q,cos(aop),sin(aop),**quiver_opts)
#fig5.savefig('AoP Cylinder, {} to {}, sun {}, quivers only.png',dpi=pps) 
#close(fig5)
