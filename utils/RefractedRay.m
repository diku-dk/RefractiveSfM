

% vi is incoming ray
% n is normal
% mu1 and mu2 are refractive indices at the boundary
% computes the outgoing refracted ray

function [vr,a,b] = RefractedRay(vi,n,mu1,mu2)

n = n/norm(n);

a = mu1/mu2;
b = -mu1*(vi'*n) - sqrt(mu1^2*(vi'*n)^2 - (mu1^2-mu2^2)*(vi'*vi));
b = b/mu2;

vr = a*vi + b*n;





