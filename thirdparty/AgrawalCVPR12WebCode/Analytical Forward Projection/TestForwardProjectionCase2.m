
% CVPR 2012 Paper Matlab Code
% Title: A Theory of Multi-Layer Flat Refractive Geometry
% Authors: Amit Agrawal, Srikumar Ramalingam, Yuichi Taguchi, Visesh Chari

% In this code, we generate simulation data for Case2. Then we use the
% forward projection equation to compute the projection of 3D points and
% show that it matches the original image projection

% Case 2 (Air-Medium-Air)

clear all;close all;clc;

mu = 1.5 ; % refractive index of medium

fprintf('Testing Forward Projection for Case 2: Air-Medium-Air\n')

fprintf('Case 2: Air-Medium-Air. Refractive index of medium = %f\n',mu);

COP = [0;0;0];  % Camera coordinate system is origin

% generate a flat refractive layer at distance d0
d0 = 15;
d1 = 24; % thickness of medium

% generate random normal of the flat refractive layer
n = [randn(1)*0.3 ; randn(1)*0.3 ; -1 + randn(1)*0.3];
n = n/sqrt(n'*n);
if(n(3)>=0)
    n = -n;
end


% display
figure;plot3(COP(1),COP(2),COP(3),'k-*');hold on;
[xx,yy] = meshgrid(-20:2:20,-20:2:20);
zz = (-d0 - n(1)*xx - n(2)*yy)/n(3);
mesh(xx,yy,zz);hold on
zz = (-(d0+d1) - n(1)*xx - n(2)*yy)/n(3);
mesh(xx,yy,zz);hold on




% Now we generate image points and find the final refracted ray working
% from the camera.
Ntrial = 100;
fprintf('Testing %d feature points\n',Ntrial);

for trial = 1:Ntrial
    
    % generate random image feature point in normalized coordinate system
    x2D = rand(2,1) - 0.5;
    
    rayD = [x2D ; 1];
    
    %find intersection of this ray with the flat refracting plane
    lamda = -d0/(n'*rayD);
    
    % lamda should be > 0 so that 3D point is in front
    while(lamda < 0)
        x2D = rand(2,1) - 0.5;
        rayD = [x2D ; 1];
        %find intersection of this ray with the flat plane
        % lamda*vi should lie on the plane
        lamda = -d0/(n'*rayD);
    end
    
    q = lamda*rayD; % q lies on flat refracting plane
    line([COP(1) ; q(1)], [COP(2) ; q(2)],[COP(3) ; q(3)] ,'LineWidth', 2 );hold on;
    
    vi = q;  % the direction for camera ray
    
    % find refracted ray
    v2 = RefractedRay(vi,n,1,mu);
    
    % find point q2 on second layer.
    q2 = q - d1*v2/(v2'*n);
    
    line([q2(1) ; q(1)], [q2(2) ; q(2)],[q2(3) ; q(3)] ,'LineWidth', 2 );hold on;
    
    % final refracted ray will have same direction as vi
    %find a 3D point on this ray at some random distance
    fac = abs(20 + 10*rand(1));
    p = q2 + fac*vi/norm(vi);
    
    
    plot3(p(1),p(2),p(3),'g-s');hold on;
    line([p(1) ; q2(1)], [p(2) ; q2(2)],[p(3) ; q2(3)] ,'LineWidth', 2 );
    
    % now we solve the forward projection equation to find the image
    % projection of this 3D point
    M = SolveForwardProjectionCase2(d1,n,mu,p);
    
    imageProj = [M(1)/M(3); M(2)/M(3)];
    
    % the image projection should match the original feature point
    ee(trial,1) = sum(abs(imageProj - x2D));
    
end

figure;plot(ee);title('error between ground truth and computed projection using AFP');


fprintf('mean error between ground truth feature point and computed projection using AFP = %f\n',mean(ee));















