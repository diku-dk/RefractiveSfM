
% Solve Forward projection equation for two layer (Air-Medium-Air)
% tau is the thickness of medium
% p is given 3D point
% n is the normal
% mu is refractive index

% M is the 3D point on the layer closest to camera where the first
% refraction happens


function [M] = SolveForwardProjectionCase2(tau,n,mu,p)

% we can choose ANY distance for d0
d0 = tau;

d2 = d0 + tau;

M = [0;0;1];


%find POR the plane of refraction
POR = cross(n,p);
POR = POR/norm(POR);

% [z1,z2] defines a coordinate system on POR
% axis is away from the camera. z1 is along the axis
z1 = -n;
z1 = z1/norm(z1);

% find z2
z2 = cross(POR,z1);
z2 = z2/norm(z2);

% find the projection of given 3D point on POR
v = p'*z1;
u = p'*z2;


s1 = (mu - 1)*(mu + 1)*(d0 - d2 + v)^2;
s2 = (-2)*d0*u*(mu - 1)*(mu + 1)*(d0 - d2 + v);
s3 = d0^2*(d0^2*mu^2 - d0^2 - 2*d0*d2*mu^2 + 2*d0*d2 + 2*d0*mu^2*v + d2^2*mu^2 - d2^2 - 2*d2*mu^2*v + mu^2*u^2 + mu^2*v^2 - u^2);
s4 = (-2)*d0^3*mu^2*u*(d0 - d2 + v);
s5 = d0^4*mu^2*u^2;


sol = roots([s1;s2;s3;s4;s5]);

% find real solutions
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
    vi = [x;d0];
    
    v2 = RefractedRay(vi,Normal,1,mu);
    
    q2 = vi + (d0-d2)*v2/(v2'*Normal);
    
    vrd = [u;v] - q2;
    
    e = abs(vrd(1)*vi(2) - vrd(2)*vi(1));
    
    if(e < 1e-4)
        M = x*z2 + d0*z1;
        return
    end
end



















