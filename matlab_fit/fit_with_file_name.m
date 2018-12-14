%% Loading figures
clear;clc
dirname='C:\Users\Guangyuan\OneDrive\文档\UCLA 课程\Computational imaging\polarization\Pol_matlab\waterlapse2\';
% cd(dirname);
files=dir([dirname,'\*.jpg']);
j=1;
k=1;
tic
f_num=length(files);
for i=1:f_num
    a=imread([dirname,files(i).name]);% read the raw data
    stack(:,:,:,k)=imresize(a,0.3);
    figure(1);imagesc(stack(:,:,:,k)); title(num2str(i));
    drawnow;
    k=k+1;
end
toc
% rgb2gray
for i=1:size(stack,4)
stack_new(:,:,i)=double(rgb2gray(stack(:,:,:,i)));
figure(1);imagesc(stack_new(:,:,i));axis image off;drawnow;title(num2str(i));
pause(0.1)
end
% Average figures 
for i=1:4
I(:,:,i)=mean(stack_new(:,:,(10*i-9):(10*i)),3);
figure(1);imagesc(I(:,:,i));colormap jet;axis image off;title(num2str(i));
pause(0.2)
end
clear stack_new stack 
%% Caculate the Aops 
s=zeros(size(I));
s(:,:,1)=mean(I,3);
s(:,:,2)=I(:,:,3)-I(:,:,1);
s(:,:,3)=I(:,:,4)-I(:,:,2);
aop=(atan(s(:,:,3)./s(:,:,2))*0.5);
figure(1);imagesc(aop);colormap jet;axis equal;axis off;
dop=abs(sqrt(s(:,:,2).^2+s(:,:,3).^2)./s(:,:,1));
figure(2);imagesc(abs(dop));colormap jet; axis image off;
xa
aa=425:775; bb=512:811;
gra=abs(gradient(abs(aop(aa,bb))));
figure(3);imagesc(abs(gra));
[cx,cy]=find(gra=)
%% main function for estimation
% a: pixel size in the unit percent of cos(ha)
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

