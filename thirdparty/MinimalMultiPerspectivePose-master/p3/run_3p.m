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


% Gen data in the world coordinates system
clear global; clc; close all;

disp('GET THE DATA...')
[adapter,R,t,l1,l2,l3] = GetData();

% Compute the minimal solutions using 3P 
disp(' ');
disp('COMPUTE MINIMAL SOLUTIONS...')
tic
[sR,sT,sdepth1,sdepth2,sdepth3] = compute_p3(adapter);
toc

% The results for all possible solutions
disp(' ');
disp('RESULTS...');

for i=1:numel(sdepth1)
    disp(['Hyp: ',num2str(i)]);
    Ri = sR(1:3,3*(i-1)+1:3*i);
    disp(['  -> eR: (',num2str(norm(Ri-R)),')']);
    disp(['  -> et: (',num2str(norm(sT(:,i)-t)),')']);
    disp(['  -> depth1: ', num2str(sdepth1(i)), ' -> edepth1: ', num2str(l1-sdepth1(i))]);
    disp(['  -> depth2: ', num2str(sdepth2(i)), ' -> edepth2: ', num2str(l2-sdepth2(i))]);
    disp(['  -> depth3: ', num2str(sdepth3(i)), ' -> edepth3: ', num2str(l3-sdepth3(i))]);
end
