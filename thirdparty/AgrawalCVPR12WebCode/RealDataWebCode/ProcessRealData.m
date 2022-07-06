
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) MERL 2012
% CVPR 2012 Paper Title: A Theory of Multi-Layer Flat Refractive Geometry
% Author: Amit Agrawal
% Code for Real Dataset shown in Figure 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Note: This code works only for non-planar calibration object


clear all;close all;clc;

addpath(genpath('./'))

fprintf('Matlab Code for CVPR 2012 Paper\n');
fprintf('Title: A Theory of Multi-Layer Flat Refractive Geometry\n');
fprintf('Author: Amit Agrawal\n');

disp('=======================')



% Refractive index of water
mu = 1.33;
fprintf('Refractive index of water used = %f\n',mu);

dirstr = ''

% trials for RANSAC in 8 point algorithm
RANSAC_trials = 200

% image captured looking through the water tank
img_test = 5772


% CameraRays are checkerboard corner points, KK is internal calibration
% matrix of camera, Pall is the 3D points in the coordinate system of
% the left-most checkerboard
load('RealData.mat')
w = Rot2RPY(Rgt);


% %choose only right and left planes
% Pall = Pall(:,1:96);
% CameraRays = CameraRays(:,1:96);

% %choose only right and Middle planes
% Pall = Pall(:,[1:48 97:end]);
% CameraRays = CameraRays(:,[1:48 97:end]);

% %choose only left and Middle planes
% Pall = Pall(:,[49:end]);
% CameraRays = CameraRays(:,[49:end]);

% Pall = Pall(:,[97:end]);
% CameraRays = CameraRays(:,[97:end]);

N = size(Pall,2)

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontSize',16);
set(gca,'FontName','Arial');
plot3(Pall(1,:),Pall(2,:),Pall(3,:),'r*');
title('3D Points in Object Coordinate System');
grid on



fprintf('Estimating rotation and translation using perspective (central) approximation\n');
% You need to have Matlab camera calibration toolbox installed
fc = [KK(1,1);KK(2,2)];
cc = KK(1:2,3);
% [omckk,T_CA,R_CA,H,x_proj,ex,JJ] = compute_extrinsic(CameraRays(1:2,:),Pall,fc,cc,zeros(5,1),0);

% find reprojection error
% ReprojErrorCA = sqrt(mean(ex(:).^2));
% fprintf('Central Approximation: Reprojection Error = %f pixels\n',ReprojErrorCA);

% w_CA = Rot2RPY(R_CA);
% fprintf('rotation [degrees]: true = [%f,%f,%f], central approx = [%f,%f,%f]\n',w(1),w(2),w(3),w_CA(1),w_CA(2),w_CA(3));
% fprintf('translation [mm]: true = [%f,%f,%f], central approx = [%f,%f,%f]\n',Tgt(1),Tgt(2),Tgt(3),T_CA(1),T_CA(2),T_CA(3));
% 
str = sprintf('%sIMG_%04d_rect.jpg',dirstr,img_test);
im = double(imread(str));
figure;imagesc(uint8(im));hold on;
plot(CameraRays(1,:)',CameraRays(2,:)','r*','MarkerSize',8);hold on;
% plot(x_proj(1,:)',x_proj(2,:)','gs');hold on;
% str = sprintf('central approximation: Reprojection error = %f pixels\n',ReprojErrorCA);
% title(str);
% legend('detected points','reprojected 3D points');

clear XYZ omckk H x ex JJ




ray2D = inv(KK)*CameraRays;



disp('=================================')
disp('=================================')
fprintf('Now running 8 pt algorithm\n');
disp('==========8 pt =======================')

%==============8 pt ================================
solVec = EightPointEstimationTwoLayerCase2KnownMu(ray2D,Pall,mu,RANSAC_trials);
ne8pt = solVec(1:3);
tau8pt = solVec(4);
R8pt = RPY2Rot(solVec(5:7));
t8pt = solVec(8:10);
w8pt = solVec(5:7);
disp('======8 pt ==========')
fprintf('rotation [degrees]: true = [%f,%f,%f], estimated = [%f,%f,%f]\n',w(1),w(2),w(3),w8pt(1),w8pt(2),w8pt(3));
fprintf('translation [mm]: true = [%f,%f,%f], estimated = [%f,%f,%f]\n',Tgt(1),Tgt(2),Tgt(3),t8pt(1),t8pt(2),t8pt(3));
fprintf('thickness of tank in mm: ground truth = %f, estimated = %f\n',260,tau8pt);

Reproj8pt = ReprojectionErrorTwoLayersCase2(ne8pt,w8pt,t8pt,tau8pt,CameraRays,Pall,KK,mu);
fprintf('After 8 pt algo: Reprojection Error = %f pixels\n',Reproj8pt);


disp('===============================')
disp('===============================')
disp('===============================')
fprintf('We got the initial estimates. Now minimizing re-projection error using analytical forward projection equations');
%non-linear optimization
x0 = [ne8pt(1:2)/abs(ne8pt(3)) ; w8pt ; t8pt ; tau8pt ; ];
options = optimset('Display','iter','TolFun',1e-6,'TolX',1e-6,'MaxIter',5000,'MaxFunEvals',10000);
x = lsqnonlin(@CostFuncTwoLayersCase2,x0,[],[],options,CameraRays,Pall,KK,mu);
ne8pt_1 = [x(1:2) ; -1];
ne8pt_1 = ne8pt_1/norm(ne8pt_1);
w8pt_1 = x(3:5);
t8pt_1 = x(6:8);
tau8pt_1 = x(9);



disp('======After Non-Linear refinement==========')
fprintf('rotation [degrees]: true = [%f,%f,%f], estimated = [%f,%f,%f]\n',w(1),w(2),w(3),x(3),x(4),x(5));
fprintf('translation [mm]: true = [%f,%f,%f], estimated = [%f,%f,%f]\n',Tgt(1),Tgt(2),Tgt(3),x(6),x(7),x(8));
fprintf('thicknes of tank in mm: ground truth = %f, estimated = %f\n',260,x(9));
Reproj8ptBA = ReprojectionErrorTwoLayersCase2(ne8pt_1,w8pt_1,t8pt_1,tau8pt_1,CameraRays,Pall,KK,mu);
fprintf('Final reprojection error = %f pixels\n',Reproj8ptBA);



PlotReprojectedPoints(im,x0,CameraRays,Pall,mu,KK);
str = sprintf('Initial Solution after 8 pt: Reprojection error = %f pixels\n',Reproj8pt);
title(str);

PlotReprojectedPoints(im,x,CameraRays,Pall,mu,KK);
str = sprintf('Final Solution after non-linear refinement: Reprojection error = %f pixels\n',Reproj8ptBA);
title(str);








