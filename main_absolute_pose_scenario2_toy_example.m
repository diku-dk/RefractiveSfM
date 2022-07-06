%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example code for paper:
%
%       Absolute and Relative Pose Estimation in Refractive Multi View
%       In Submission to ICCV 2021
%
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the version 3 of the GNU General Public License
% as published by the Free Software Foundation.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.       
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
% IMPORTANT INFO:
% This script generates synthetic data for absolute pose estimation under
% refraction. 
%
clc;clear all;close all;

addpath utils/
addpath solvers\
addpath libso3\

K = [4800 0 960;0 4800 540;0 0 1];

stdnoise = 1;
opt_z = [rand(1,2)*0.5-0.25,-1]';opt_z = opt_z ./ norm(opt_z);
rotaxis = cross([0;0;1],opt_z);
rotaxis = rotaxis ./ norm(rotaxis);
theta = acos(dot(opt_z,[0;0;1]));
R = expSO3(theta*rotaxis);
assert(abs(dot(R*[0;0;1],opt_z)-1)<1e-6);
% camera pose
R1 = R; C1 = [rand(1,2),rand(1)+1]';
refra_n1 = 1.0;% air refractive index;
refra_n2 = 1.5;% water refractive index;
miu = refra_n2/refra_n1;%

n1w = [0;0.5;1];n1w = n1w./norm(n1w);
d1w = -0.5;
num_features = 100;

isplanar = 0;

[q,Q,ref_As]=ab_dataGenerator(num_features,R1,C1,n1w,d1w,refra_n1,refra_n2,stdnoise,K,isplanar);

figure;colormap jet;
plot3(C1(1),C1(2),C1(3),'k-d','MarkerSize',15);hold on;
plot3(0,0,0,'m-s','MarkerSize',15);hold on;
[xx,yy] = meshgrid(min(ref_As(1,:))-1:0.1:max(ref_As(1,:))+1,min(ref_As(2,:))-1:0.1:max(ref_As(2,:))+1);
zz = (-d1w - n1w(1)*xx - n1w(2)*yy)/n1w(3);
h3=mesh(xx,yy,zz);hold on

h1=plotCameraFrustum(R1,-R1*C1,K,0.5,'r',1.5);
for ii = 1:num_features
    h=line([C1(1) ; ref_As(1,ii)], [C1(2) ; ref_As(2,ii)],[C1(3) ; ref_As(3,ii)] ,'LineWidth', 2 , 'Color', 'b');hold on;
    h.Color(4) = 0.5;
    h2=plot3(ref_As(1,ii),ref_As(2,ii),ref_As(3,ii),'y-s','MarkerFaceColor','y');hold on;
    h=line([Q(1,ii) ; ref_As(1,ii)], [Q(2,ii) ; ref_As(2,ii)],[Q(3,ii) ; ref_As(3,ii)] ,'LineWidth', 2, 'Color', 'g');
    h.Color(4) = 0.5;
end
h4=plot3(Q(1,:),Q(2,:),Q(3,:),'ko','MarkerFaceColor','k');hold on;
xlabel('x: (m)','Interpreter','latex');ylabel('y: (m)','Interpreter','latex');zlabel('z: (m)','Interpreter','latex');
title('Simulation case: regular','Interpreter','latex');
legend([h1,h2,h3,h4],{'Camera','Refractive Point','Refractive Plane','World Point'},'Interpreter','latex');

%% Proposed solution.
addpath(genpath('./thirdparty/AgrawalCVPR12WebCode/RealDataWebCode'))   
tic;
[RPO,CPO] = PnP_Refraction_scenario2(Q,q,K,n1w,d1w,refra_n1,refra_n2,isplanar);
duration = toc;
rmpath(genpath('./thirdparty/AgrawalCVPR12WebCode/RealDataWebCode'))    

disp(['PO runtime:', num2str(duration)]);
disp(['PO Rotation error:', num2str(norm(RPO'*R1-eye(3)))]);
disp(['PO Location error:', num2str(norm(CPO-C1))]);

%% Argawal linear solution.
tic;
addpath(genpath('./thirdparty/AgrawalCVPR12WebCode/RealDataWebCode'))   
[Ragw,Cagw] = abs_pose_agw(Q,q,K,refra_n1,refra_n2,0,isplanar,0);
rmpath(genpath('./thirdparty/AgrawalCVPR12WebCode/RealDataWebCode'))    
duration = toc;
disp(['AGWL runtime:', num2str(duration)]);
disp(['AGWL Rotation error:', num2str(norm(Ragw'*R1-eye(3)))]);
disp(['AGWL Location error:', num2str(norm(Cagw-C1))]);

%% Argawal Iterative solution.
tic;
addpath(genpath('./thirdparty/AgrawalCVPR12WebCode/RealDataWebCode'))   
[RagwI,CagwI] = abs_pose_agw(Q,q,K,refra_n1,refra_n2,1,isplanar,0);
rmpath(genpath('./thirdparty/AgrawalCVPR12WebCode/RealDataWebCode'))    
duration = toc;
disp(['AGW runtime:', num2str(duration)]);
disp(['AGW Rotation error:', num2str(norm(RagwI'*R1-eye(3)))]);
disp(['AGW Location error:', num2str(norm(CagwI-C1))]);

%% Hanner 5pt
addpath('./thirdparty/refractive_pose-master')
tic
[Rana,Cana] = Hanner_8pt(Q,q,K,miu,R1,C1);
duration = toc;
disp(['5pt runtime:', num2str(duration)]);
disp(['5pt Rotation error:', num2str(norm(Rana'*R1-eye(3)))]);
disp(['5pt Location error:', num2str(norm(Cana-C1))]);
rmpath('./thirdparty/refractive_pose-master')

%% EPnP: depends on mexopencv
try
    tic
    [rvec,tvec] = cv.solvePnP(Q',q',K,'Method','EPnP');
    duration = toc;
    REPnP = expSO3(rvec);
    CEPnP = -REPnP'*tvec;
    disp(['EPnP runtime:', num2str(duration)]);
    disp(['EPnP Rotation error:', num2str(norm(REPnP'*R1-eye(3)))]);
    disp(['EPnP Location error:', num2str(norm(CEPnP-C1))]);
catch
    REPnP = eye(3);
    CEPnP = zeros(3,1);
    warning('Failure: You need to install mexopencv in order to run this section.');
end

figure
cmap = jet(5);
plot3(0,0,0,'r.');
hold on;

h1=plotCameraFrustum(R1',-R1'*C1,K,2,cmap(1,:),2);
h2=plotCameraFrustum(RPO',-RPO'*CPO,K,2,cmap(2,:),2);
h3=plotCameraFrustum(Ragw',-Ragw'*Cagw,K,2,cmap(3,:),2);
h4=plotCameraFrustum(RagwI',-RagwI'*CagwI,K,2,cmap(4,:),2);
h5=plotCameraFrustum(REPnP',-REPnP'*CEPnP,K,2,cmap(5,:),2);

view(3);
axis equal;
title('Relative Pose', 'FontName', 'Arial', 'FontSize', 10);
legend([h1,h2,h3,h4,h5],{'Truth','Proposed','AGW','AGWI','5pt'}, 'FontName', 'Arial', 'FontSize', 10);
grid minor
