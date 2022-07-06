
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) MERL 2012
% CVPR 2012 Paper Title: A Theory of Multi-Layer Flat Refractive Geometry
% Author: Amit Agrawal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [] = PlotReprojectedPoints(im,x,CameraRays,Pall,mu,KK)

N = size(Pall,2);

n = [x(1:2) ; -1];
n = n/norm(n);


figure;imagesc(uint8(im));hold on;
plot(CameraRays(1,:)',CameraRays(2,:)','r*');hold on;
R1 = RPY2Rot(x(3:5));
t1 = x(6:8);
tau = x(9);
P = R1*Pall + t1*ones(1,N);

clear xreproj

% find the image reprojection point for each 3D point
for ii = 1:N
    p = P(:,ii);
    % solve analytical forward projection equation
    M = SolveFlatRefractionPointTwoLayerCase2(tau,n,mu,p);
    M = KK*M;
    xreproj(1,ii) = M(1)/M(3);
    xreproj(2,ii) = M(2)/M(3);
end
plot(xreproj(1,:)',xreproj(2,:)','gs');hold on;

legend('detected points','reprojected 3D points');


