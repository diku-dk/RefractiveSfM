function varargout = forward_project(n,d,n1,n2,t,M)
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    refra_n = n2 / n1;
    % M = [2.045114831690593,0,1]';
    O = t + d*n;%[t2(1),t2(2),0]';
    Mproj = t + dot(M-t,n)*n;
    x = norm(Mproj-M);
    z = norm(Mproj-O);

    %% solve quartic 
    u = solve_u(1/refra_n,d,z,x);
 
    r1 = M - Mproj;
    r1 = r1 ./ norm(r1);
    Q = O + r1.*u;

    i1 = Q - M; 
    i1 = i1./norm(i1); % incident ray 1
    nn = -n;
    miu = n1/n2;%
    tr1 = miu * i1 + sqrt(1-miu^2*(1-(dot(nn,i1))^2)) * nn - miu * dot(nn,i1) * nn;% transmitted ray 1
    q1 = -tr1 ./ norm(tr1);
    
    varargout{1} = q1;
    varargout{2} = i1;
    varargout{3} = Q;
end

function u = solve_u(n,d,z,x)
% b = [ u^4, u^3, u^2, u, 1]
    a = [ n^2 - 1, -2*x*(n^2 - 1), x^2*(n^2 - 1) - z^2 + d^2*n^2, -2*d^2*n^2*x, d^2*n^2*x^2];
    u = roots(a);
    valid = abs(imag(u)) < 1e-6;
    u = real(u(valid));
    valid = u > 0 & u < x;
    u = u(valid);
end