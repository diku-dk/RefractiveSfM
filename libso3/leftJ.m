
function J = leftJ(r)
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
    tolerance = 1e-12;
    angle = norm(r);
    if angle < tolerance
        % If the angle is small, fall back on the series representation
        N = 10;
        J = eye(3);
        pxn = eye(3);
        px = hat(r);
        for n = 1:N
            pxn = pxn*px/(n + 1);    
            J = J + pxn;
        end
    else
        axis = r/angle;

        cph = (1 - cos(angle))/angle;
        sph = sin(angle)/angle;

        J = sph * eye(3) + (1 - sph) * axis * axis' + cph * hat(axis);
    end       
end