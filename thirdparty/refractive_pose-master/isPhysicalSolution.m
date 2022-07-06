%[b reason] = isPhysicalSolution(P,X,x,N,W,r) returns true if
%	camera matrix P back-projects image points x through the refractive
%	plane passing through the point W with normal N such that the rays
%	intersect scene points X. r is the refractive index ratio n1/n2.
%
%	'reason' contains any failure code:
%	1: Camera center and scene point on same side of plane.
%	2: Image ray does not intersect the plane.
%	3: Total internal reflection occurred at the interface.
%	4: The refracted ray and scene point are on opposite sides of the plane.
%	5: The refracted ray did not pass within the given tolerance of the
%	   scene point.


function [b reason] = isPhysicalSolution(P,X,x,N,W,r,tol)

if nargin<7
	tol = 1e-4;
end
b = false;
reason = 0;

N = N/norm(N);
C = -P(:,1:3)'*P(:,4);
u = P(:,1:3)'*[x(1:2,:); ones(1,size(x,2))];
u = bsxfun(@rdivide,u,sqrt(sum(u.^2)));
u2 = zeros(size(u));
P2 = zeros(size(u));

for i=1:size(u,2)
	if sign((C-W)'*N)==sign((X(:,i)-W)'*N), reason = 1; return; end
	ct1 = -(N'*u(:,i));
	t = N'*(C-W)/ct1; if t<0, reason = 2; return; end
	P2(:,i) = C + t*u(:,i);
	ct2 = 1 - r^2*(1-ct1^2); if ct2<0, reason = 3; return; end
	ct2 = sign(ct2)*sqrt(abs(ct2));
	u2(:,i) = r*u(:,i) + (r*ct1-sign(ct1)*ct2)*N; 
end

for i=1:size(u,2)
	v = X(:,i) - P2(:,i);
	if v'*u2(:,i)<0, reason = 4; return; end
	d = norm(v-(v'*u2(:,i))*u2(:,i));
	if d>tol, reason = 5; return; end
end
b = true;
