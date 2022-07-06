

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) MERL 2012
% CVPR 2012 Paper Title: A Theory of Multi-Layer Flat Refractive Geometry
% Author: Amit Agrawal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the total reprojection error for all 3D points using estimated
% calibration parameteres

function [ee] = ReprojectionErrorTwoLayersCase2(n,w,t,tau,rayDnoise,XYZ,KK,mu)

n = n/norm(n);
R = RPY2Rot(w);


N = size(XYZ,2);
C = zeros(2*N,1);


for ii = 1:N
    
    p = R*XYZ(:,ii) + t;
    
    M = SolveFlatRefractionPointTwoLayerCase2(tau,n,mu,p);
    M = KK*M;

    px1 = M(1)/M(3);
    py1 = M(2)/M(3);
    
    px = rayDnoise(1,ii);
    py = rayDnoise(2,ii);
    
    C(2*ii-1,1) = px-px1;
    C(2*ii,1) = py-py1;
    
end

C = C(:);


ee = sqrt(mean(C.^2));






