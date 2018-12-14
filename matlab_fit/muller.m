function [MM]=muller(theta_i,n)
theta_r=asin(sin(theta_i)/n);
e1=cos(theta_i-theta_r)^2+1;
e2=cos(theta_i-theta_r)^2-1;
e3=2*cos(theta_i-theta_r);
M=[e1,e2,0,0;
    e2,e1,0,0;
    0,0,e3,0;
    0,0,0,e3];
pre=(sin(2*theta_i)/sin(2*theta_r))/(2*sin(theta_i+theta_r)^2/...
    cos(theta_i-theta_r)^2);
MM=pre*M;
end