function [adapter,R,C,l1,l2]=GetData(Plot,Noise_Level)
% AUTORIGHTS
% ---------------------------------------------------------
% Copyright (c) 2018, Pedro Miraldo
% 
% This file is part of the Minimal Multi-Perspective Pose
% code and is available under the terms of the MIT License
% provided in LICENSE. Please retain this notice and
% LICENSE if you use this file (or any portion of it)
% in your project.
% ---------------------------------------------------------


if nargin == 0
    Plot = 1;
    Noise_Level = 0;
end

K = [9.842439e+02, 0.000000e+00, 6.900000e+02;
     0.000000e+00, 9.808141e+02, 2.331966e+02;
     0.000000e+00, 0.000000e+00, 1.000000e+00];

dangle = 2;

 
% get rand data
C = [randn()*5-10; randn()*5; randn()*5];
a1 = randn()*pi/9; a2 = randn()*pi/9; a3 = randn()*pi/9;
R = angle2dcm(a1,a2,a3);

P1  = [randn()*2.5-15;randn()*5+15;randn()*5];
P2  = [randn()*2.5-10;randn()*5-15;randn()*5];
LP1W = [randn()*2.5;randn()*2.5-7.5;randn()*2.5];
LP2W = [randn()*2.5;randn()*2.5+7.5;randn()*2.5];
L1  = [LP1W,LP2W];

C1C = [0;5;0]+randn(3,1)*2;
C1W = R*C1C + C;
D1W = P1-C1W; D1W = D1W/norm(D1W);
D1C = R'*D1W; D1C = D1C/norm(D1C);
Rcamera1 = [1,0,0;0,0,-1;0,1,0]*angle2dcm(randn()*dangle*pi/180,randn()*dangle*pi/180,randn()*dangle*pi/180);
if Noise_Level ~= 0
    u   = K*Rcamera1*D1C;
    u   = u/u(3);
    D1C = K\(u + [randn(2,1)*Noise_Level; 0]);
    D1C = Rcamera1'*(D1C/norm(D1C));
end

C2C = [2;-2;2]+randn(3,1)*2;
C2W = R*C2C + C;
D2W = P2-C2W; D2W = D2W/norm(D2W);
D2C = R'*D2W; D2C = D2C/norm(D2C);
Rcamera2 = [1,0,0;0,0,1;0,-1,0]*angle2dcm(randn()*dangle*pi/180,randn()*dangle*pi/180,randn()*dangle*pi/180);
if Noise_Level ~= 0
    u   = K*Rcamera2*D2C;
    u   = u/u(3);
    D2C = K\(u + [randn(2,1)*Noise_Level; 0]);
    D2C = Rcamera2'*(D2C/norm(D2C));
end

LP1C = R'*LP1W -R'*C;
LP2C = R'*LP2W -R'*C;
LP1C = LP1C/norm(LP1C);
LP2C = LP2C/norm(LP2C);
Rcamera0 = [1,0,0;0,1,0;0,0,1]*angle2dcm(randn()*dangle*pi/180,randn()*dangle*pi/180,randn()*dangle*pi/180);
if Noise_Level ~= 0
    u   = K*Rcamera2*LP2C;
    u   = u/u(3);
    LP2C = K\(u + [randn(2,1)*Noise_Level/2; 0]);
    LP2C = Rcamera2'*LP2C/norm(LP2C);
    
    u   = K*Rcamera2*LP1C;
    u   = u/u(3);
    LP1C = K\(u + [randn(2,1)*Noise_Level/2; 0]);
    LP1C = Rcamera2'*(LP1C/norm(LP1C));
end
N1C = cross(LP1C,LP2C);
N1C = N1C/norm(N1C);

% output data
l1 = norm(P1-C1W);
l2 = norm(P2-C2W);

% set the adapter
% camera info
adapter.camera.c1 = C1C;
adapter.camera.c2 = C2C;
adapter.camera.d1 = D1C;
adapter.camera.d2 = D2C;
adapter.camera.l1.normal = N1C;
% world info
adapter.world.p1 = P1;
adapter.world.p2 = P2;
adapter.world.l1.p1 = LP1W;
adapter.world.l1.p2 = LP2W;

% plot data
if(Plot ~= 0)
    
    close all
    
    figure(1);
    hold on;
    quiver3(0,0,0,3,0,0,'Color',[0.6350, 0.0780, 0.1840],'LineWidth',2);
    quiver3(0,0,0,0,3,0,'Color',[0.4660, 0.6740, 0.1880],'LineWidth',2);
    quiver3(0,0,0,0,0,3,'Color',[0, 0.4470, 0.7410],'LineWidth',2);
    text(0,0,-1,'W');
    
    plot3(L1(1,:),L1(2,:),L1(3,:),'--b','LineWidth',2);
    plot3(C(1),C(2),C(3),'ok','LineWidth',1,'MarkerSize',5);
    text(C(1),C(2),C(3)-1,'C')
    fill3([L1(1,1),L1(1,2),C(1)],[L1(2,1),L1(2,2),C(2)],[L1(3,1),L1(3,2),C(3)], 'r', 'FaceAlpha', 0.5);
    ec1 = R*[1,0,0]';
    ec2 = R*[0,1,0]';
    ec3 = R*[0,0,1]';
    quiver3(C(1),C(2),C(3),3*ec1(1),3*ec1(2),3*ec1(3),'Color',[0.6350, 0.0780, 0.1840],'LineWidth',2);
    quiver3(C(1),C(2),C(3),3*ec2(1),3*ec2(2),3*ec2(3),'Color',[0.4660, 0.6740, 0.1880],'LineWidth',2);
    quiver3(C(1),C(2),C(3),3*ec3(1),3*ec3(2),3*ec3(3),'Color',[0, 0.4470, 0.7410],'LineWidth',2);
    
    plot3(P1(1),P1(2),P1(3),'ob','LineWidth',1,'MarkerSize',5);
    text(P1(1),P1(2),P1(3)-1,'P1')
    plot3(C1W(1),C1W(2),C1W(3),'or','LineWidth',1,'MarkerSize',5);
    text(C1W(1),C1W(2),C1W(3)-1,'C1');
    
    quiver3(C1W(1),C1W(2),C1W(3),3*D1W(1),3*D1W(2),3*D1W(3),'Color','r','LineWidth',2);
    plot3([C1W(1),P1(1)],[C1W(2),P1(2)],[C1W(3),P1(3)],':r');
    
    plot3(P2(1),P2(2),P2(3),'ob','LineWidth',1,'MarkerSize',5);
    text(P2(1),P2(2),P2(3)-1,'P2')
    plot3(C2W(1),C2W(2),C2W(3),'or','LineWidth',1,'MarkerSize',5);
    text(C2W(1),C2W(2),C2W(3)-1,'C2');
    quiver3(C2W(1),C2W(2),C2W(3),3*D2W(1),3*D2W(2),3*D2W(3),'Color','r','LineWidth',2);
    plot3([C2W(1),P2(1)],[C2W(2),P2(2)],[C2W(3),P2(3)],':r');
    
    MID_PL = (LP1W+LP2W)/2;
    DIR_L  = (LP2W-LP1W); DIR_L = DIR_L/norm(DIR_L);
    plot3([MID_PL(1)+20*DIR_L(1),MID_PL(1)-20*DIR_L(1)],...
        [MID_PL(2)+20*DIR_L(2),MID_PL(2)-20*DIR_L(2)],...
        [MID_PL(3)+20*DIR_L(3),MID_PL(3)-20*DIR_L(3)],...
        '--b','LineWidth',2);
    
    
    view(-120,35);
    grid on;
    box on;
    axis equal;
    hold off;
    
end

end
