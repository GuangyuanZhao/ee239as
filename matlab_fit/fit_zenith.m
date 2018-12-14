    
%% main function for estimation
% a: pixel size in the unit percent of cos(ha)
% dir='C:\Users\Guangyuan\OneDrive\文档\UCLA 课程\Computational imaging\polarization\Pol_matlab\data_important';
% cd(dir);
clear ;
load('dop_new.mat');  load('aop_new.mat'); load('aop_horizon.mat');
maximum=2;
a=maximum/301;

sizeyx = size (aop_new);
gra=abs(gradient(aop_new));
[cx,cy]=find(dop_new==min(dop_new(:)));
cx=146; % 主要是因为dop 不够准，我目测了一个
% [x, y] = meshgrid (1:sizeyx (2), 1:sizeyx (1));
[x, y] = meshgrid (1:sizeyx (2), 1:sizeyx (1));

% [x, y] = meshgrid (0.5:sizeyx (2)-.5, .5:sizeyx (1)-.5);
x = (x-cx)*a;
y = (y-cy)*a;
[theta,r]=cart2pol(x,y);
%% Do the fitting using optimization 
% Initial value for the fitting
phi_s=deg2rad(133);
h_s=deg2rad(30);
a=35/301*2/301; %认为301 个是相当于r=1,我们实际只有35个 pixel 
params0=[phi_s,h_s];
plot_on=0;
b=0.1;
cx=146;cy=142;
c=[cx,cy];
data=abs(aop_new);
[paramsF] = fit_aop (data, params0, a, b,c, plot_on); % Use a simpled model for the fitting. 

elevation_fit=rad2deg(paramsF(2));
azimuth_fit=rad2deg(paramsF(1));

