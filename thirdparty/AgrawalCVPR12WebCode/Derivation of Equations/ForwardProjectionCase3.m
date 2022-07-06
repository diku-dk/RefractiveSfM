

% CVPR 2012 Paper Matlab Code
% Title: A Theory of Multi-Layer Flat Refractive Geometry
% Authors: Amit Agrawal, Srikumar Ramalingam, Yuichi Taguchi, Visesh Chari


% This code computes the forward projection equation for two layer system
% (Air-Medium1-Medium2) (Case 3). The resulting equation is 12th degree

% Camera is assumed to be in air with Refractive index = 1
% Refractive index of medium = mu1
% Refractive index of last medium = mu2
% distance of first layer from camera = d0
% Thickness of medium1 = d1

% The analysis is done on the plane of refraction.


clear all;close all;

disp('Forward Projection Equation for Case 3')

syms d0 d1 mu1 mu2 real  % known calibration parameters
syms ux uy real          % known projection of 3D point P on the plane of refraction
syms x real              % The unknown camera ray can be parameterized as [x;d_0]

u = [ux ; uy];

% normal (axis) in the 2D plane
n = [0;-1];

vp0 = [x;d0];  % unknown camera ray

c0 = vp0.'*n;

%==========================================
% find refracted vector vp1 = a1*vp0 + b1*n
a1 = 1/mu1;
b1 = -c0 - sqrt(c0^2 - (1-mu1^2)*(vp0.'*vp0));
b1 = b1/mu1;

disp('Finding refracted vector vp1');
vp1 = a1*vp0 + b1*n;

disp('Finding refracted vector vp2');
% since (vp0.'*vp0) = (vp1.'*vp1), we will use (vp0.'*vp0) instead of
% (vp1.'*vp1)

a2 = mu1/mu2;
b2 = -mu1*(vp1.'*n) - sqrt(mu1^2*(vp1.'*n)^2 - (mu1^2-mu2^2)*(vp0.'*vp0));
b2 = b2/mu2;
vp2 = a2*vp1 + b2*n;

q1 = vp0;
q2 = q1 - d1*vp1/(vp1.'*n);


% Forward projection equation
% take cross product of vp2 and u-q2

disp('Forward Projection Equation: cross(vp2,u-q2)=0')

t = u - q2;

FP = vp2(2)*t(1) - vp2(1)*t(2);
[nn,~] = numden(FP);
FP = simplify(nn);


%remove the square root term
AA = c0^2 - (1-mu1^2)*(vp0.'*vp0);
AA = expand(AA);

BB = mu1^2*(vp1.'*n)^2 - (mu1^2-mu2^2)*(vp0.'*vp0);
BB = expand(BB);


ca = coeffs(FP,sqrt(AA));
cb = coeffs(ca(1),sqrt(BB));
cc = coeffs(ca(2),sqrt(BB));

k1 = cc(1)
k2 = cc(2)
k3 = cb(1)

%check
%simplify(k1*sqrt(AA) + k2*sqrt(AA)*sqrt(BB) + k3*sqrt(BB) - FP)
%equation is given by  k1*sqrt(AA) + k2*sqrt(AA)*sqrt(BB) + k3*sqrt(BB) = 0
% solution after remvoing sqrt terms is

FP = (k1^2*AA + k3^2*BB +  - k2^2*AA*BB)^2 - 4*k1^2*k3^2*AA*BB;

disp('forward projection equation')
pretty(FP)

disp('degree of forward projection equation for x')
[tt1,tt2] = coeffs(FP,x);
tt2







