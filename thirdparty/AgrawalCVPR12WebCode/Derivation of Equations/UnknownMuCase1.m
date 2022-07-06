

% CVPR 2012 Paper Matlab Code
% Title: A Theory of Multi-Layer Flat Refractive Geometry
% Authors: Amit Agrawal, Srikumar Ramalingam, Yuichi Taguchi, Visesh Chari


% This code derives the 6th degree equation to compute the calibration
% parameters for Case 1 (Air-Medium) under unknown refractive index (after
% axis estimation).
% We work on the plane of refraction itself

clear all;close all;clc

disp('Deriving 6th degree equation for Case 1')

syms d0 mu1 real  % unknown calibration parameters
syms alpha real   % unknown translation magnitude along the axis
syms ux uy real   % known projection of 3D point R*P + t_A^{\perp} on the plane of refraction
syms vx vy real   % known camera ray vp0 = [vx ; vy]
syms gama real

u = [ux ; uy + alpha];

% normal (axis) in the 2D plane
n = [0;-1];

vp0 = [vx;vy];  % known camera ray. Since camera ray is known, we can assume it is normalized, so vp0.*vp0= 1

c0 = vp0.'*n;

%==========================================
% find refracted vector vp1 = a1*vp0 + b1*n
a1 = 1/mu1;
b1 = -c0 - sqrt(c0^2 - (1-mu1^2));
b1 = b1/mu1;

disp('Finding refracted vector vp1');
vp1 = a1*vp0 + b1*n;

q1 = -d0*vp0/c0;



% take cross product of vp1 and u - q1
t = u - q1;

eq = vp1(2)*t(1) - vp1(1)*t(2);
[nn,~] = numden(eq);
eq = simplify(nn);

%remove the square root term
AA = c0^2 - (1-mu1^2);
AA = expand(AA);

c = coeffs(eq,sqrt(AA));
eq = c(2)^2*AA - c(1)^2;
eq = simplify(eq);


disp('FRC equation')
pretty(eq)

disp('substituting gama as mu1^2')


eq = subs(eq,mu1,sqrt(gama));



disp('FRC equation')
pretty(eq)


disp('generating two more equations for two more correspondences')
syms vx2 vy2 ux2 uy2 real
syms vx3 vy3 ux3 uy3 real

eq2 = subs(eq,vx,vx2);
eq2 = subs(eq2,vy,vy2);
eq2 = subs(eq2,ux,ux2);
eq2 = subs(eq2,uy,uy2);
pretty(eq2)

eq3 = subs(eq,vx,vx3);
eq3 = subs(eq3,vy,vy3);
eq3 = subs(eq3,ux,ux3);
eq3 = subs(eq3,uy,uy3);
pretty(eq3)


disp('solving for gama')
gsol = solve(eq,gama);
gsol = simplify(gsol);

disp('substituting gama in Eq2 and Eq3')

eq2 = subs(eq2,gama,gsol);
[nn,~] = numden(eq2);
eq2 = simplify(nn);
clear nn

eq3 = subs(eq3,gama,gsol);
[nn,~] = numden(eq3);
eq3 = simplify(nn);
clear nn


disp('powers of alpha and d0 in eq2')
[tt1,tt2] = coeffs(eq2,d0);tt2
[tt1,tt2] = coeffs(eq2,alpha);tt2


disp('removing alpha^2 from eq2 and eq3')
c2 = coeffs(eq2,alpha);
c3 = coeffs(eq3,alpha);


eqAlpha = c2(3)*eq3 - c3(3)*eq2;
eqAlpha = factor(eqAlpha);
eqAlpha = simplify(eqAlpha/(d0*vx-vy*ux)^2);

disp('obtained linear equation in alpha')

%now eqAlpha is a linear equation in alpha
% find coeffs

k = coeffs(eqAlpha,alpha);

%eqAlpha = k(2) * alpha + k(1). Thus, alpha can be solved as -k(1)/k(2)

% substitue back into eq3. eq3 can be written as
%eq3 = c3(3)*alpha^2 + c3(2)*alpha + c3(1)

disp('obtained alpha')

disp('substiuting alpha into eq3')
f2 = c3(3)*k(1)^2 - c3(2)*k(1)*k(2) + c3(1)*k(2)^2;


% this gives a 8th degree equation which can be factoriezed into a 6th
% degree equation and two linear terms given by
% (d0*vx*vx3*vy - d0*vx*vx3*vy3 + ux*vx3*vy*vy3 - ux3*vx*vy*vy3)
% (d0*vx*vx3*vy + d0*vx*vx3*vy3 - ux*vx3*vy*vy3 - ux3*vx*vy*vy3)

ff = (d0*vx*vx3*vy - d0*vx*vx3*vy3 + ux*vx3*vy*vy3 - ux3*vx*vy*vy3)*(d0*vx*vx3*vy + d0*vx*vx3*vy3 - ux*vx3*vy*vy3 - ux3*vx*vy*vy3);

disp('Simplifying: This could take up to 15-20 mins depending on your machine')

tic
f2 = factor(f2);
finaleq = simplify(f2/ff);
toc


disp('degree of equation for d0')
[tt1,tt2] = coeffs(finaleq,d0);tt2





























