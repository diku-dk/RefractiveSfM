
function r = logSO3(R)
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
    r = [0;0;0];
    res = (trace(R)-1)*0.5;
    if res < -1
        res = -1;
    end
    angle = acos(res);
    if angle > 1e-10
        so3 = angle / (2*sin(angle))*(R-R');
        r = [-so3(2,3);so3(1,3);-so3(1,2)];
    else
        so3 = zeros(3,3);
        for i = 1:2
            so3 = so3 + (-1)^(i-1)/i.*(R-eye(3))^(i);
        end
        so3 = 0.5.*(so3-so3');
        r = [-so3(2,3);so3(1,3);-so3(1,2)];
    end
end