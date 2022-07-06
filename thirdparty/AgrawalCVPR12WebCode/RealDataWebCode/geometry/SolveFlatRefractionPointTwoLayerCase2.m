
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) MERL 2012
% CVPR 2012 Paper Title: A Theory of Multi-Layer Flat Refractive Geometry
% Author: Amit Agrawal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% solve the forward projection for flat refraction and compute the point on
% the water surface where refraction happens

% equation of water surface is n'*x + d = 0;
% P is the 3D point in the camera coordinate system
% this is actually the way to compute the reprojection point of a 3D point. need to redrive to single plane case.

function [M] = SolveFlatRefractionPointTwoLayerCase2(tau,n,mu,p)

% we can choose any value of d!
d = tau;

d2 = d + tau;

M = zeros(3,1);

%refraction plane
RefPlane = cross(n,p);
RefPlane = RefPlane/norm(RefPlane);

% take first axis to be normal direction, away from the camera
z1 = -n;
z1 = z1/norm(z1);

z2 = cross(RefPlane,z1);
z2 = z2/norm(z2);

v = p'*z1;
u = p'*z2;


%solve 4th degree equation
s1 = (mu - 1)*(mu + 1)*(d - d2 + v)^2;
s2 = (-2)*d*u*(mu - 1)*(mu + 1)*(d - d2 + v);
s3 = d^2*(d^2*mu^2 - d^2 - 2*d*d2*mu^2 + 2*d*d2 + 2*d*mu^2*v + d2^2*mu^2 - d2^2 - 2*d2*mu^2*v + mu^2*u^2 + mu^2*v^2 - u^2);
s4 = (-2)*d^3*mu^2*u*(d - d2 + v);
s5 = d^4*mu^2*u^2;


sol = roots([s1;s2;s3;s4;s5]);

idx = find(abs(imag(sol)) < 1e-6);
if(isempty(idx))
    disp('no solution');
    return
end

sol1 = sol(idx);
nn = size(sol1,1);


Normal = [0;-1];

for ii = 1:nn

    x = sol1(ii,1);
    vi = [x;d];

    v2 = RefractedRay(vi,Normal,1,mu);

    q2 = vi + (d-d2)*v2/(v2'*Normal);



    vrd = [u;v] - q2;


    e = abs(vrd(1)*vi(2) - vrd(2)*vi(1));

    if(e < 1e-4)
        M = x*z2 + d*z1;
        return
    end
end



















