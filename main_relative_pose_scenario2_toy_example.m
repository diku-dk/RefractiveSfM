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

K = [4800 0 960;0 4800 540;0 0 1];

stdnoise = 0.5;
% camera pose
refra_n1 = 1.0;% refractive index;
refra_n2 = 1.5;% refractive index;
miu = refra_n2/refra_n1;%

n1w = [0;0.5;1];n1w = n1w./norm(n1w);
d1w = -0.5;
num_features = 100;

R1 = eye(3);% from c to w
n1w = [rand(2,1)*0.5-0.25;1];n1w = n1w./norm(n1w);
d1w = -(2 + rand(1));

COP1 = [0,0,0]';
COP2 = [rand(1,2),rand(1)-3]' + [2,0,0]';

% 0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0
redo = 0;% sometimes, the assertion will fail due to the very bad geometry
while true
%     try
        [imgf1s,imgf2s,ref_As,ref_Bs,Qs,R2] = rel_dataGenerator(num_features,R1,n1w,d1w,COP1,COP2,refra_n2,refra_n1,stdnoise,K);
%     catch
%         redo = 1;
%     end
    if redo == 0, break; end
end

%% display
cmap = lines(9);
figure;
plot3(COP1(1),COP1(2),COP1(3),'-d','Color',cmap(1,:),'MarkerSize',15);hold on;
h1=plotCameraFrustumRC(R1',COP1,K,1,cmap(1,:),1.5);
plot3(COP2(1),COP2(2),COP2(3),'-d','Color',cmap(2,:),'MarkerSize',15);hold on;
h2=plotCameraFrustumRC(R2',COP2,K,1,cmap(2,:),1.5);
[xx,yy] = meshgrid(min([ref_As(1,:), ref_Bs(1,:)])-1:0.1:max([ref_As(1,:), ref_Bs(1,:)])+1,...
                   min([ref_As(2,:), ref_Bs(2,:)])-1:0.1:max([ref_As(2,:), ref_Bs(2,:)])+1);
zz = (-d1w - n1w(1)*xx - n1w(2)*yy)/n1w(3);
h3=mesh(xx,yy,zz);hold on

for ii = 1:num_features
    line([COP1(1) ; ref_As(1,ii)], [COP1(2) ; ref_As(2,ii)],[COP1(3) ; ref_As(3,ii)] ,'LineWidth', 2 , 'Color', cmap(3,:));hold on;
    h5=plot3(ref_As(1,ii),ref_As(2,ii),ref_As(3,ii),'-o','Color',cmap(4,:),'MarkerFaceColor',cmap(4,:));hold on;
    line([Qs(1,ii) ; ref_As(1,ii)], [Qs(2,ii) ; ref_As(2,ii)],[Qs(3,ii) ; ref_As(3,ii)] ,'LineWidth', 2, 'Color', cmap(5,:));

    line([COP2(1) ; ref_Bs(1,ii)], [COP2(2) ; ref_Bs(2,ii)],[COP2(3) ; ref_Bs(3,ii)] ,'LineWidth', 2, 'Color', cmap(6,:));hold on;
    h6=plot3(ref_Bs(1,ii),ref_Bs(2,ii),ref_Bs(3,ii),'-o','Color',cmap(7,:),'MarkerFaceColor',cmap(7,:));hold on;
    line([Qs(1,ii) ; ref_Bs(1,ii)], [Qs(2,ii) ; ref_Bs(2,ii)],[Qs(3,ii) ; ref_Bs(3,ii)] ,'LineWidth', 2, 'Color', cmap(8,:));
end
h4=plot3(Qs(1,:),Qs(2,:),Qs(3,:),'ko','MarkerFaceColor','k');hold on;
xlabel('x: (m)','Interpreter','latex');ylabel('y: (m)','Interpreter','latex');zlabel('z: (m)','Interpreter','latex');
title('Simulation case: relative pose','Interpreter','latex');
legend([h1,h2,h3,h4,h5,h6],{'Camera 1','Camera 2','Refractive Plane','World Point','Refractive Point 1','Refractive Point 2'},'Interpreter','latex');
% end

Rt = R1'*R2;
Ct = COP2;
Xtrue = Qs;

%% Nister five point
try
tic
 E1 = cv.findEssentialMat(imgf1s(1:2,:)', imgf2s(1:2,:)', 'CameraMatrix', K);
[Rord, tord, good, mask, xyzPoints] = cv.recoverPose(E1, imgf1s(1:2,:)', imgf2s(1:2,:)', ...
    'CameraMatrix', K, 'DistanceThreshold', 10);
R5pt = Rord';
C5pt = -Rord'*tord;
duration = toc;
disp(['Nister five point runtime:', num2str(duration)]);
disp(['Nister five point Rotation error:', num2str(norm(R5pt'*Rt-eye(3)))]);
catch
    warning('Failure: You need to install mexopencv in order to run this section.');
end

%% proposed: initialization using Nister five point
try
tic;
[X_ve,R_ve,C_ve,time_ve] = refinement_virtual_epipolar_scenario2(n1w, d1w, refra_n1, refra_n2, imgf1s, imgf2s, K);
duration = toc;
% recorver real trans
disp(['proposed (Nister five point) runtime:', num2str(duration)]);
disp(['proposed (Nister five point) Rotation error:', num2str(norm(R_ve'*Rt-eye(3)))]);
catch
    warning('Failure: You need to install mexopencv in order to run this section.');
end

%% proposed: initialization using truth + perturbation, no dependency on mexopencv
std = 0.1;
R0 = Rt * expSO3(randn(3,1)*std);
C0 = Ct + randn(3,1)*std;
tic;
[X_ve,R_ve,C_ve,time_ve] = refinement_virtual_epipolar_scenario2(n1w, d1w, refra_n1, refra_n2, imgf1s, imgf2s, K, R0, C0);
duration = toc;
% recorver real trans
disp(['proposed (perturbation) runtime:', num2str(duration)]);
disp(['proposed (perturbation) Rotation error:', num2str(norm(R_ve'*Rt-eye(3)))]);