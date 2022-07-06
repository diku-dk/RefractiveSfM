

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) MERL 2012
% CVPR 2012 Paper Title: A Theory of Multi-Layer Flat Refractive Geometry
% Author: Amit Agrawal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% estimate tz and depth for two layers given known mu, axis, R and n
% cross t
% case 1 mu2 = mu0
% do two point ransac

function [d,t] = EstimateTzTwoLayeCase2KnownMuRansac(n,R,s,mu,xy,XYZ,idx)


N = size(XYZ,2);


% find rotation to make axis along z axis
w = cross(n,[0;0;-1]);
w = w/norm(w);
theta = acos(n'*[0;0;-1]);
w = theta*w;
Rc = w2R(w);
clear w
clear n




ss = Rc * s;
% now tx and tz can be obtained. We need to find tz
t = [-ss(2) ; ss(1) ; 0];
clear ss


% rotate 3D points
XYZ = Rc*R*XYZ + t*ones(1,N);

n2D = [0;-1];

% find z1, z2
z1 = [0;0;1];

%
zp = [0;1];
if(isempty(idx))
    keyboard
end


Nsmall = size(idx,1);
A = zeros(Nsmall,2);
b = zeros(Nsmall,1);

for jj = 1:Nsmall
    
    ii = idx(jj,1);
    
    % camera ray
    v = Rc*xy(:,ii);
    
    P = XYZ(:,ii);
    
    %find z2
    z2 = cross(z1,cross(z1,v));
    z2 = z2/norm(z2);
    % make along x axis
    if(z2(1)<0)
        z2 = -z2;
    end
    
    %find projection
    vp = [z2'*v ; z1' * v];
    vp = vp/norm(vp);
    
    %find refracted ray
    [v1,~,~] = RefractedRay(vp,n2D,1,mu);
    
    
    % final refracted ray is vp itself
    u = [z2'*P ; z1'*P];
    
    %make linear system
    b(jj,1) = -Cross2(vp,u);
    
    A(jj,:) = [Cross2(vp,v1)/(v1'*n2D) Cross2(vp,zp)];
    
end


idx1 = nchoosek([1:Nsmall],2);

nn = size(idx1,1);

bestEr = Inf;
for ii = 1:nn
    vecid = idx1(ii,:);
    vecid = vecid(:);
    x = A(vecid,:)\b(vecid);
    d = x(1);
    tz = x(2);
    er = median((A*x - b).^2);
    if(er < bestEr)
        bestEr = er;
        bestDTz = [d;tz];
    end
end

d = bestDTz(1);
tz = bestDTz(2);


% put tz in t
t(3) = tz;
% rotate back into original coordinate system
t = Rc'*t;
















