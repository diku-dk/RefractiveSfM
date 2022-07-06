
function vechat = hat(vec)
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
    if size(vec,1) == 3 
        vechat = [  0,     -vec(3),  vec(2);
                vec(3),   0    , -vec(1);
               -vec(2),  vec(1),   0    ];  
    elseif size(vec,1) == 6
        vechat = [ hat( vec(4:6,1) ) vec(1:3,1); zeros(1,4) ];    
    end   
end