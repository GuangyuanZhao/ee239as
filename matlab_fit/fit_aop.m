function [paramsF,E] = fit_aop (data, params0, a, b,c,plot_on)
% Input Parameters:
% Required:
% data -- the aop image.
% params0 -- The user's initial guess of the parameters to be fit:
%           [phi_s,h_s]
%           phi_s,h_s should be specified in radians.
%           note: as the fit minimizes only localy it is important for this
%           inital guess to be fairly close to the true value.
% a -- pixel size in unit cos(h_p)       
% b -- border if we want to compare in a circular region around the center
% c -- center of the aop map
% Optional:
% plot_on -- a binary parameter that determines whether the outcome is
% plotted.

sizeyx = size (data);
[x, y] = meshgrid (0.5:sizeyx (2)-.5, .5:sizeyx (1)-.5);
x = (x-c(1))*a;
y = (y-c(2))*a;
[theta,r]=cart2pol(x,y);


% The funtion to be minimized is the negative of the log likelihood
datafun = @(params)(sum (sum ((expected (theta,r,params,b))))...
                        -sum (sum (data.*log (expected (theta,r,params,b)))));
options = optimset ('MaxFunEvals', 10000, 'MaxIter', 100, 'TolFun', 1e-2);
% fminsearch performs the multivariable minimization
[paramsF,fval,exitflag,output]  = fminsearch (datafun, params0, options);
E=(expected (theta,r,paramsF,b));
end

function E = expected (theta,r,params,b)
% (x,y,params,a)
% The expected aop per pixel.
[m, n]=size(r);
alpha=zeros(m,n);
phi_s=params(1);
h_s=params(2);

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
            alpha(i,j)=abs(atan(tana));
        end
    end
end


b;
E=alpha;
end



