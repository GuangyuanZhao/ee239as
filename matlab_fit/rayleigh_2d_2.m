function [x,y,aop]=rayleigh_2d_2(Hs,As)
	
% 	@Hs: elevation angle of the sun (0,90)
% 	@As: azimuth angle of the sun (-180,180)

	% maximum dgree of polarization in Rayleigh model
	deltamax = 1;
	% normalized radius of celestial sphere
	R = 1;


% 	# thetas: zenith angle of the sun (90-Hs)
	thetas = pi/2 - Hs*pi/180.0;
	As = As*pi/180;
% 	# phis: azimuth angle of the sun (0,360)
% 	# relationship between sun and the sphere coordinate
	if pi/2 <= As <= pi
		phis = As - pi/2;
    else
		phis = As + 1.5*pi;
    end

% 	# dop: degree of polarization
% 	# aop: angle of polarization	 
	dop = zeros(100,400) ;% # np.zeros(), dtype=float;;
	aop = zeros(100,400);
	x = zeros(100,400);
	y =zeros(100,400);

% 	# initial location of the sun
	xs=0;ys=0;
% 	# for loop	 
	i=1;j=1;
Ho=linspace(0,pi/2,100);
Ao=linspace(0,2*pi,400);
% 	# Ho: elevation angle of the observer
% 	# Ao: azimuth angle of the observer
     for i=1:length(Ho)
        tic
		for j=1:length(Ao)
			thetao = pi/2 - Ho(i);
			phio = Ao(j) ;% # for simplification;
% 			# if pi/2 <= Ao <= pi:
% 			# 	phio = Ao - pi/2
% 			# else:
% 			# 	phio = Ao + 1.5*pi
% 			# gamma: angular difference between sky point and the sun
			gamma = acos(sin(thetas)*sin(thetao)*cos(phio-phis)+cos(thetas)*cos(thetao));

% 			# dop: degree of polarization
			dop(i,j)=(deltamax*(sin(gamma))^2)/(1+cos(gamma)^2);

% 			# aop: angle of polarization
			if sin(phio-phis)*sin(thetas) == 0
				aop(i,j) = pi/2;
            else
				aop(i,j) = atan((sin(thetao)*cos(thetas)-cos(thetao)*cos(phio-phis)*sin(thetas))/(sin(phio-phis)*sin(thetas)));
            end
			x(i,j) = R*sin(thetao)*cos(phio);
			y(i,j) = R*sin(thetao)*sin(phio);
        end
        toc
    end
	aop =rad2deg(aop);

% 	# location of the sun
% 	# directly sphere mapping
	xs = R*sin(thetas)*cos(phis);
	ys = R*sin(thetas)*sin(phis);


	meridian_x = zeros(2);
	meridian_y = zeros(2);

	meridian_x(1) =xs/(sqrt(xs^2 + ys^2));
	meridian_x(2)=meridian_x(2)- meridian_x(1);

	meridian_y(1)= ys/(sqrt(xs^2 + ys^2));
	meridian_y(2) =meridian_y(2) -meridian_y(1);

end