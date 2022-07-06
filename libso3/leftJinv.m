
function Jinv = leftJinv(r)
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
    phi = r;
    ph = norm(phi);
    if ph < tolerance
        % If the angle is small, fall back on the series representation        
        Jinv = eye(3);
        pxn = eye(3);
        px = hat(phi);
        for n = 1:10
            pxn = pxn * px/n;
            Jinv = Jinv + lut_bernoullinumber(n) * pxn;
        end
    else
        axis = phi/norm(phi);
        ph_2 = 0.5*ph;

        Jinv =   ph_2 * cot(ph_2)* eye(3)...
               + (1 - ph_2 * cot(ph_2))* axis * axis'...
               - ph_2 * hat(axis);
    end   
end