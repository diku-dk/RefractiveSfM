function [sR,sT] = compute_3l(adapter)
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



% Get data from the adapter
% Camera data
C1           = adapter.Camera.C1;
C2           = adapter.Camera.C2;
IterPlane_L0 = adapter.Camera.IterPlane_L0;
IterPlane_L1 = adapter.Camera.IterPlane_L1;
IterPlane_L2 = adapter.Camera.IterPlane_L2;
% world data
P01 = adapter.World.L0.P1;
P02 = adapter.World.L0.P2;
P11 = adapter.World.L1.P1;
P12 = adapter.World.L1.P2;
P21 = adapter.World.L2.P1;
P22 = adapter.World.L2.P2;

% compute the null space
uc1 = IterPlane_L0(1:3)/norm(IterPlane_L0(1:3));
uc2 = IterPlane_L1(1:3)/norm(IterPlane_L1(1:3));
uc3 = IterPlane_L2(1:3)/norm(IterPlane_L2(1:3));

vw1 = (P02-P01); vw1 = vw1/norm(vw1);
vw2 = (P12-P11); vw2 = vw2/norm(vw2);
vw3 = (P22-P21); vw3 = vw3/norm(vw3);
uw1 = cross(P01,vw1);
uw2 = cross(P11,vw2);
uw3 = cross(P21,vw3);

uc = [uc1,uc2,uc3];
vw = [vw1,vw2,vw3];
uw = [uw1,uw2,uw3];

C = [zeros(3,1), C1, C2];

B_ = [ uc1(1)*vw1(1), uc1(2)*vw1(1), uc1(3)*vw1(1), uc1(1)*vw1(2), uc1(2)*vw1(2), uc1(3)*vw1(2), uc1(1)*vw1(3), uc1(2)*vw1(3), uc1(3)*vw1(3);
       uc2(1)*vw2(1), uc2(2)*vw2(1), uc2(3)*vw2(1), uc2(1)*vw2(2), uc2(2)*vw2(2), uc2(3)*vw2(2), uc2(1)*vw2(3), uc2(2)*vw2(3), uc2(3)*vw2(3);
       uc3(1)*vw3(1), uc3(2)*vw3(1), uc3(3)*vw3(1), uc3(1)*vw3(2), uc3(2)*vw3(2), uc3(3)*vw3(2), uc3(1)*vw3(3), uc3(2)*vw3(3), uc3(3)*vw3(3)];
[~,~,V] = svd(B_);
B = V(:,4:9);

[b1,b2,b3,b4,b5] = solvecoeffs(B(:,1),B(:,2),B(:,3),B(:,4),B(:,5),B(:,6));

% compute the results
nsolns = length(b1);
sR = zeros(3,nsolns);
sT= zeros(3,nsolns);

A = zeros(9,4);

for i = 1:nsolns
    
    % get the rotation parameters
    b = [b1(i);b2(i);b3(i);b4(i);b5(i);1];
    P = B*b;
    myr=P(1:9);
    myR=reshape(myr,3,3);
    a = norm(myR(1,:));
    if det(myR)<0; a=-a; end
    myR = (myR')/a;

    % estimate the translation for the respective rotation
    r1 = myR(1,1); r2 = myR(1,2); r3 = myR(1,3);
    r4 = myR(2,1); r5 = myR(2,2); r6 = myR(2,3);
    r7 = myR(3,1); r8 = myR(3,2); r9 = myR(3,3);
    
    for ii = 1 : 3
        
        c1 = C(1,ii); c2 = C(2,ii); c3 = C(3,ii);
        uc1 = uc(1,ii); uc2 = uc(2,ii); uc3 = uc(3,ii);
        vw1 = vw(1,ii); vw2 = vw(2,ii); vw3 = vw(3,ii);
        uw1 = uw(1,ii); uw2 = uw(2,ii); uw3 = uw(3,ii);
        
        A((ii-1)*3+1:ii*3,:) = [vw1*(r2*uc2 + r3*uc3) + vw2*(r5*uc2 + r6*uc3) + vw3*(r8*uc2 + r9*uc3),                                       -uc2*(r1*vw1 + r4*vw2 + r7*vw3), -uc3*(r1*vw1 + r4*vw2 + r7*vw3), - vw1*(uc2*(c1*r2 - c2*r1) + uc3*(c1*r3 - c3*r1)) - vw2*(uc2*(c1*r5 - c2*r4) + uc3*(c1*r6 - c3*r4)) - vw3*(uc2*(c1*r8 - c2*r7) + uc3*(c1*r9 - c3*r7)) - uw1*(r2*uc3 - r3*uc2) - uw2*(r5*uc3 - r6*uc2) - uw3*(r8*uc3 - r9*uc2);
            -uc1*(r2*vw1 + r5*vw2 + r8*vw3), vw1*(r1*uc1 + r3*uc3) + vw2*(r4*uc1 + r6*uc3) + vw3*(r7*uc1 + r9*uc3), -uc3*(r2*vw1 + r5*vw2 + r8*vw3),   vw1*(uc1*(c1*r2 - c2*r1) - uc3*(c2*r3 - c3*r2)) + vw2*(uc1*(c1*r5 - c2*r4) - uc3*(c2*r6 - c3*r5)) + vw3*(uc1*(c1*r8 - c2*r7) - uc3*(c2*r9 - c3*r8)) + uw1*(r1*uc3 - r3*uc1) + uw2*(r4*uc3 - r6*uc1) + uw3*(r7*uc3 - r9*uc1);
            -uc1*(r2*vw1 + r5*vw2 + r8*vw3), vw1*(r1*uc1 + r3*uc3) + vw2*(r4*uc1 + r6*uc3) + vw3*(r7*uc1 + r9*uc3), -uc3*(r2*vw1 + r5*vw2 + r8*vw3),   vw1*(uc1*(c1*r2 - c2*r1) - uc3*(c2*r3 - c3*r2)) + vw2*(uc1*(c1*r5 - c2*r4) - uc3*(c2*r6 - c3*r5)) + vw3*(uc1*(c1*r8 - c2*r7) - uc3*(c2*r9 - c3*r8)) + uw1*(r1*uc3 - r3*uc1) + uw2*(r4*uc3 - r6*uc1) + uw3*(r7*uc3 - r9*uc1)];
        
    end
    As = A(:,1:3);
    bs = A(:,4);
    test = -inv(As'*As)*As'*bs; 
    
    % store the results
    sR(:,(i-1)*3+1:i*3) = myR;
    sT(:,i) = -myR*test;
    
end
