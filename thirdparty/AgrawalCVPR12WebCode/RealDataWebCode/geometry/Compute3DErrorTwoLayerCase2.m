

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) MERL 2012
% CVPR 2012 Paper Title: A Theory of Multi-Layer Flat Refractive Geometry
% Author: Amit Agrawal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute 3D error between the back-projected ray and the given 3D point
% for two layers Case 2
% done, note that this is quite special case

function [C] = Compute3DErrorTwoLayerCase2(t,xy,XYZ,R,n,mu,tau)


% assume some d1
d1 = tau;

N = size(XYZ,2);

% transform points
P = R*XYZ + t*ones(1,N);

e = zeros(N,1);

for ii = 1:N

    pt = P(:,ii);
    v = xy(:,ii);

    q1 = -d1*v/(v'*n);


    [v1,~,~] = RefractedRay(v,n,1,mu);


    %find q2

    q2 = q1 - tau*v1/(v1'*n);


    %final refracted ray is v, this is special case for case 2 where miu0 = miu2, so after two refraction, v2 = v

    e(ii,1) = norm(cross(v,q2 - pt))/norm(v);

end


C = median(e.^2);

