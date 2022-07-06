
function vpt = para_trans(R1,R2,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% support code for paper:
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
%     A = logSO3(R1'*R2);
%     Rpt = expSO3(A/2);
%     vpt = invhat(R2'*R1*Rpt*hat(v)*Rpt);
    % fast version
    vpt = expSO3(logSO3(R2'*R1)/2)*v;
end


%function vpt = para_trans(R1,R2,v)
%    A = logSO3(R1'*R2);
%    Rpt = expSO3(A/2);
%    if size(v,2) == 1
%        vpt = invhat(R2'*R1*Rpt*hat(v)*Rpt);
%    else
%        vpt = invhat(R2'*R1*Rpt*(R1'*v)*Rpt);
%    end
%end