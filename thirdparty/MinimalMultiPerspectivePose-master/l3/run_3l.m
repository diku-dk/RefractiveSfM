% AUTORIGHTS
% ---------------------------------------------------------
% Copyright (c) 2018, Pedro Miraldo
% 
% This file is part of the Minimal Multi-Perspective Pose
% code and is available under the terms of the MIT License
% provided in LICENSE. Please retain this notice and
% LICENSE if you use  this file (or any portion of it)
% in your project.
% ---------------------------------------------------------


% Gen data in the world coordinates system
clear global; clc; close all;

disp('GET THE DATA...')
[R,C,adapter] = GetData();

% Compute the minimal solutions using 3L
disp(' ');
disp('COMPUTE MINIMAL SOLUTIONS...')
tic
[sR,sT] = compute_3l(adapter);
toc

% The results for all possible solutions
disp(' ');
disp('RESULTS...');

for i=1:numel(sT(1,:))
    disp(['Hyp: ',num2str(i)]);
    Ri = sR(1:3,3*(i-1)+1:3*i);
    disp(['  -> eR: (',num2str(norm(Ri-R)),')']);
    disp(['  -> et: (',num2str(norm(sT(:,i)-C)),')']);
end
