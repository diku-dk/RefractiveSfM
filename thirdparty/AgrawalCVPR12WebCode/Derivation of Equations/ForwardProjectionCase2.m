

% CVPR 2012 Paper Matlab Code
% Title: A Theory of Multi-Layer Flat Refractive Geometry
% Authors: Amit Agrawal, Srikumar Ramalingam, Yuichi Taguchi, Visesh Chari


% This code computes the forward projection equation for two layer system
% (Air-Medium-Air) (Case 2). The resulting equation is 4th degree

% Camera is assumed to be in air with Refractive index = 1
% Refractive index of medium = mu1
% distance of first layer from camera = d0
% Thickness of medium = d1

% The analysis is done on the plane of refraction.


clear all;close all;clc;

disp('Forward Projection Equation for Case 2')

syms d0 d1 mu1 real  % known calibration parameters
syms ux uy real      % known projection of 3D point P on the plane of refraction
syms x real          % The unknown camera ray can be parameterized as [x;d_0]

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


q1 = vp0;

q2 = q1 - d1*vp1/(vp1.'*n);

% Forward projection equation
% take cross product of vp0 and u-q2
% since final refracted ray is parallel to vp0

disp('Forward Projection Equation: cross(vp0,u-q2)=0')

t = u - q2;

FP = vp0(2)*t(1) - vp0(1)*t(2);
[nn,~] = numden(FP);
FP = simplify(nn);


%remove the square root term
AA = c0^2 - (1-mu1^2)*(vp0.'*vp0);
AA = expand(AA);

c = coeffs(FP,sqrt(AA));
FP = c(2)^2*AA - c(1)^2;
FP = simplify(FP);

disp('forward projection equation')
pretty(FP)

disp('degree of forward projection equation for x')
[tt1,tt2] = coeffs(FP,x);
tt2







