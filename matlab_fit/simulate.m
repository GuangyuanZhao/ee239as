%% Simulate the the skylight aop above or beneath the water
%------------------------------------------------
% Aop : angle of polarization
%% Built the coordinates
maximum=2; %2r
px=301; % Sampling of simulation 
a=maximum/px;
x=1:px; y=1:py;
cx=(px+1)/2; cy=(py+1)/2;
x = (x-cx)*a;
y = (y-cy)*a;
[theta,r]=cart2pol(x,y); 
%%
%-----------------------------------------------
% Input sun's position, in degree 
phi_s=deg2rad(90);
h_s=deg2rad(30);
%-----------------------------------------------

[m, n]=size(r);
alpha=zeros(m,n);
for i=1:m
    for j=1:n
        phi_p=theta(i,j);
        cosh_p=r(i,j);
        sinh_p=sqrt(1-cosh_p^2);
        % Calculate the distribution of tan(a) using equation 
        if sin(phi_p-phi_s)*cos(h_s) == 0
			alpha(i,j) = pi/2;
        else
            tana=(cosh_p*sin(h_s)-sinh_p*cos(h_s)*cos(phi_s-phi_p))/(cos(h_s)*sin(phi_s-phi_p));
            alpha(i,j)=atan(tana);
        end
    end
end

% generate a mask for simulate the sun 
rr=cos(h_s);
index_r=find((r<rr+0.04)&(r>rr-0.1));
index_p=find((theta>phi_s-0.04)&(theta<phi_s+0.1));
[c]=intersect(index_r,index_p);
mask2=ones(size(alpha));
mask2(c)=0;

% generate a mask 
mask=(r<=1);
figure(1);imagesc(abs(mask_2.*mask.*alpha));colormap jet; axis equal;axis off;
figure(2);imagesc((aop_new));colormap jet;axis equal; axis off;

%% Calculate Dop &I
Dop=zeros(m,n);
I=zeros(m,n);

Dop_max=0.95;
I_max=1;
thetas=pi/2-h_s;
phis=phi_s;
for i=1:m
    for j=1:n
        phi_p=theta(i,j);
        h_p=acos(r(i,j));
        thetao=pi/2-h_p;
        gamma = acos(sin(thetas)*sin(thetao)*cos(phi_p-phis)+cos(thetas)*cos(thetao));
        % Dop: degree of polarization
        Dop(i,j) = (Dop_max*(sin(gamma))^2)/(1+cos(gamma)^2);
        % Intensity
        I(i,j)=I_max/2*(1+cos(gamma)^2);
    end
end
figure(4);imagesc(abs(Dop));
%% The whole stokes of the sky 

S1=I;
S2=I.*Dop.*cos(2.*(alpha));
S3=I.*Dop.*sin(2.*(alpha));
S4=zeros(size(S1));

S2(c)=0; % corresponds to the S2 of sun
S3(c)=0; % corresponds to the S3 of sun

%Encapuslate the components into the s. I Q U V 
s=cat(3,S1,S2,S3,S4);
figure(5);imagesc(abs(I));

%% water rafraction using Fresnel refraction function 
s_1=zeros(size(s));
nindex=1.33;
Dop_final=zeros(size(Dop));
Aop_final=zeros(size(alpha)); % Aop after refraction 

for i=1:m
    for j=1:n
        phi_p=theta(i,j);
        thetao=asin(r(i,j));
        s_1(i,j,:)=(muller(abs(thetao),nindex))*reshape(s(i,j,:),[4,1]); % Fresnel refraction 
        Dop_final(i,j)=sqrt(s_1(i,j,2)^2+s_1(i,j,3)^2)/s_1(i,j,1);
        Aop_final(i,j)=0.5*atan2(real(s_1(i,j,3)),real(s_1(i,j,2)));
    end
end

%% Display 
figure(5);imagesc(abs(Dop_final)); colormap jet;
figure(6);imagesc(abs(Aop_final));colormap jet;
figure(7);imagesc(abs(mask.*alpha));colormap jet;
figure(8);imagesc(Aop_final-alpha);colormap jet;  % Difference map 