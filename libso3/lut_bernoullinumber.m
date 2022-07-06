
function b = lut_bernoullinumber(n)
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
    lut = [1,-0.500000000000000,0.166666666666667,0,-0.0333333333333333,0, ...
           0.0238095238095238,0,-0.0333333333333334,0,0.0757575757575749,0, ...
           -0.253113553113554,0,1.16666666666667,0,-7.09215686274515,0,54.9711779448621,0,-529.124242424242];
    b = lut(n+1);
end