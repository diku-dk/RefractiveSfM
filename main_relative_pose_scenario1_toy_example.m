%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright 2022, Xiao Hu, Francois Lauze, Kim Steenstrup Pedersen
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License. 
%   You may obtain a copy of the License at
%
%       http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.
%
% IMPORTANT INFO:
% This script generates synthetic data for relative pose estimation under
% refraction. 
%
clc;clear all;close all;

addpath utils/
addpath solvers\
addpath libso3\

NN = 100;
aa = 8*pi/64;
R1 = [cos(aa) 0 sin(aa);0 1 0;-sin(aa) 0 cos(aa)];% from c to w
aa = pi/16;
R2 = [cos(-aa) 0  sin(-aa);0 1 0;-sin(-aa) 0 cos(-aa)];% from c to w
ior1 = 1;
ior2 = 1.5;
n1w = R1 * [0,0,1]';d1 = 0;% assume normal of refractive plane is parallel to the optical axis
n2w = R2 * [0,0,1]';d2 = 0;% assume normal of refractive plane is parallel to the optical axis
COP1 = R1*[0,0,-2]' + [0,0,0]';
COP2 = R2*[0,0,-3]' + [2,0,0]';
d2 = -3 - n2w'*COP2;

KK = [4800           0         960; ...
         0        4800         540; ...
         0           0           1];

noise_level = 0;
[imgf1s,imgf2s,ref_As,ref_Bs,Qs] = rel_scenario1_dataGenerator(NN,R1',R2',n1w,n2w,COP1,COP2,ior2,ior1,d1,d2,noise_level,KK);

     
%% display
cmap = lines(9);
figure;
plot3(COP1(1),COP1(2),COP1(3),'-d','Color',cmap(1,:),'MarkerSize',15);hold on;
h1=plotCameraFrustum(R1',COP1,KK,1,cmap(1,:),1.5);
plot3(COP2(1),COP2(2),COP2(3),'-d','Color',cmap(2,:),'MarkerSize',15);hold on;
h2=plotCameraFrustum(R2',COP2,KK,1,cmap(2,:),1.5);
[xx,yy] = meshgrid(min(ref_As(1,:))-1:0.1:max(ref_As(1,:))+1,min(ref_As(2,:))-1:0.1:max(ref_As(1,:))+1);
zz = (-d1 - n1w(1)*xx - n1w(2)*yy)/n1w(3);
h3=mesh(xx,yy,zz);hold on
[xx,yy] = meshgrid(min(ref_Bs(1,:))-1:0.1:max(ref_Bs(1,:))+1,min(ref_Bs(2,:))-1:0.1:max(ref_Bs(1,:))+1);
zz = (-d2 - n2w(1)*(xx) - n2w(2)*yy)/n2w(3);
mesh(xx,yy,zz);hold on
for i = 1:NN
    line([COP1(1) ; ref_As(1,i)], [COP1(2) ; ref_As(2,i)],[COP1(3) ; ref_As(3,i)] ,'LineWidth', 2 , 'Color', cmap(3,:));hold on;
    h5=plot3(ref_As(1,i),ref_As(2,i),ref_As(3,i),'-o','Color',cmap(4,:),'MarkerFaceColor',cmap(4,:));hold on;
    line([Qs(1,i) ; ref_As(1,i)], [Qs(2,i) ; ref_As(2,i)],[Qs(3,i) ; ref_As(3,i)] ,'LineWidth', 2, 'Color', cmap(5,:));
    
    line([COP2(1) ; ref_Bs(1,i)], [COP2(2) ; ref_Bs(2,i)],[COP2(3) ; ref_Bs(3,i)] ,'LineWidth', 2, 'Color', cmap(6,:));hold on;
    h6=plot3(ref_Bs(1,i),ref_Bs(2,i),ref_Bs(3,i),'-o','Color',cmap(7,:),'MarkerFaceColor',cmap(7,:));hold on;
    line([Qs(1,i) ; ref_Bs(1,i)], [Qs(2,i) ; ref_Bs(2,i)],[Qs(3,i) ; ref_Bs(3,i)] ,'LineWidth', 2, 'Color', cmap(8,:));
end
h4=plot3(Qs(1,:),Qs(2,:),Qs(3,:),'ko','MarkerFaceColor','k');hold on;
xlabel('x: (m)','Interpreter','latex');ylabel('y: (m)','Interpreter','latex');zlabel('z: (m)','Interpreter','latex');
title('Simulation case: relative pose','Interpreter','latex');
legend([h1,h2,h3,h4,h5,h6],{'Camera 1','Camera 2','Refractive Plane','World Point','Refractive Point 1','Refractive Point 2'},'Interpreter','latex');
%% end

%% new section
lambda = ior1/ior2;
miu = ior2/ior1;
d1c = abs(dot(n1w,COP1)+d1);
d2c = abs(dot(n2w,COP2)+d2);

%%%%%%%%%%%%%%%%%%%%%% begin check optimization performance %%%%%%%%%%%%%%%
% do some virtual camera transformation
n1c = [0,0,1]';
n2c = [0,0,1]';

%% ground truth
t1 = -R1'*COP1;% t in c
t2 = -R2'*COP2;% t in c
R = R1'*R2; % 2->1
C = -R1'*R2*t2 + t1;% 2->1
Xtrue = R1'*(Qs-COP1);

%% first, use pless solver to find out the relative pose
% pless solver only works for 0 - 0.1, from 0.2, the rotation is very bad
tic
[R17pt,t17pt] = pless_solver_opengv_17pt(imgf1s,imgf2s,KK,n1c,n2c,d1c,d2c,ior1,ior2,R);% providing R only for selecting the best sols from multiple candidates
toc

tic
[Rlin,tlin] = refractive_rel_pose_linear_solver(imgf1s,imgf2s,KK,n1c,n2c,d1c,d2c,ior1,ior2,R);% providing R only for selecting the best sols from multiple candidates
toc

%% pinhole ignoring refraction.
E1 = cv.findEssentialMat(imgf1s(1:2,:)', imgf2s(1:2,:)', 'CameraMatrix', KK);
[Rord, tord, good, mask, xyzPoints] = cv.recoverPose(E1, imgf1s(1:2,:)', imgf2s(1:2,:)', ...
    'CameraMatrix', KK, 'DistanceThreshold', 10);

%% how to initialize the structure
R0 = Rord';
C0 = -Rord'*tord;

figure;
hold on;
h1=plotCameraFrustum(R',C,KK,2,cmap(2,:),2);
h2=plotCameraFrustum(R17pt',C,KK,1,cmap(3,:),1.5);
h3=plotCameraFrustum(Rlin',C,KK,1,cmap(4,:),1.5);
h4=plotCameraFrustum(R0',C,KK,1,cmap(5,:),1.5);
view(3);
axis equal;
title('Relative Pose Initialization', 'FontName', 'Arial', 'FontSize', 10);
legend([h1,h2,h3,h4],{'Cam2: Truth','Cam2: LHD','Cam2: Proposed','Cam2: Essential Decomposition'}, 'FontName', 'Arial', 'FontSize', 10);
grid minor

%%%%%%%%%%%%%%%%%%%%%% begin check optimization performance %%%%%%%%%%%%%%%
% chadebecq iccv using virtual camera 
% [X_cbq_iccv_vc,R_cbq_iccv_vc,C_cbq_iccv_vc,timelsq_vc] = refinement_chadebecq_iccv_vc(imgf1s,imgf2s,n1c,n2c,d1c,d2c,ior1,ior2,miu,KK,R0,C0);
% s_cbq_iccv_vc = mean(vecnorm(diff(X_cbq_iccv_vc')') ./ vecnorm(diff(Xtrue')'));
% C_cbq_iccv_vc = C_cbq_iccv_vc./s_cbq_iccv_vc;

% iccv using forward projection
[X_cbq_iccv,R_cbq_iccv,C_cbq_iccv,timelsq] = refinement_chadebecq_iccv(imgf1s,imgf2s,n1c,n2c,d1c,d2c,ior1,ior2,miu,KK,R0,C0);
s_cbq_iccv = mean(vecnorm(diff(X_cbq_iccv')') ./ vecnorm(diff(Xtrue')'));
C_cbq_iccv = C_cbq_iccv./s_cbq_iccv;


% chadebecq ijcv using virtual camera 
% [X_cbq_ijcv_vc,R_cbq_ijcv_vc,C_cbq_ijcv_vc,time_cbq_ijcv_vc] = refinement_chadebecq_ijcv_vc(imgf1s,imgf2s,n1c,n2c,d1c,d2c,ior1,ior2,miu,KK,R0,C0);
% C_cbq_ijcv_vc = (C_cbq_ijcv_vc + [0;0;d1c] - R_cbq_ijcv * [0;0;d2c]);
% s_cbq_ijcv_vc = mean(vecnorm(diff(X_cbq_ijcv_vc')') ./ vecnorm(diff(Xtrue')'));
% C_cbq_ijcv_vc = C_cbq_ijcv_vc./s_cbq_ijcv_vc;

% ijcv using forward projection
[X_cbq_ijcv,R_cbq_ijcv,C_cbq_ijcv,time_cbq_ijcv] = refinement_chadebecq_ijcv(imgf1s,imgf2s,n1c,n2c,d1c,d2c,ior1,ior2,miu,KK,R0,C0);
% recorver real trans
C_cbq_ijcv = (C_cbq_ijcv + [0;0;d1c] - R_cbq_ijcv * [0;0;d2c]);
s_cbq_ijcv = mean(vecnorm(diff(X_cbq_ijcv')') ./ vecnorm(diff(Xtrue')'));
C_cbq_ijcv = C_cbq_ijcv./s_cbq_ijcv;

% vec
[X_ve,R_ve,C_ve,time_ve] = refinement_virtual_epipolar(imgf1s,imgf2s,n1c,n2c,d1c,d2c,ior1,ior2,KK,R0,C0);
% recorver real trans
C_ve = (C_ve + [0;0;d1c] - R_ve * [0;0;d2c]);
s_ve = mean(vecnorm(diff(X_ve')') ./ vecnorm(diff(Xtrue')'));
C_ve = C_ve./s_ve;


figure
plot3(Xtrue(1,:),Xtrue(2,:),Xtrue(3,:),'ro');hold on;grid minor;

X_cbq_iccv = X_cbq_iccv./s_cbq_iccv;
plot3(X_cbq_iccv(1,:),X_cbq_iccv(2,:),X_cbq_iccv(3,:),'b*');hold on;grid minor;
% % 
X_cbq_ijcv = X_cbq_ijcv./s_cbq_ijcv;
plot3(X_cbq_ijcv(1,:),X_cbq_ijcv(2,:),X_cbq_ijcv(3,:),'cs');hold on;grid minor;

X_ve = X_ve./s_ve;
plot3(X_ve(1,:),X_ve(2,:),X_ve(3,:),'gd');hold on;grid minor;
 
legend({'Ground Truth','Chad-ICCV','Chad-IJCV','Optimized-Virt'});



