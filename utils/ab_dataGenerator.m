function varargout = ab_dataGenerator(N,R1,C1,nw,dw,refra_n1,refra_n2,varargin)
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
% Inputs:
%   N: number of features.
%   R1: True camera rotation from the camera frame to world frame.
%   C1: Camera location in thr world frame.
%   nw,dw: parameters of refractive interface.
%   refra_n1,refra_n2: refractive indices.
% Optional:
%   std: Noise standard deviation.
%   K: Camera Intrinsic Matrix.
%   isplanar: 1 or 0.
%
% Outputs:
%   q1s: feature points.
%   Q1s: 3D points.
%   ref_As: refractive points.
%
    if nargin == 7
        std = 0*1e-6;
    else
        std = varargin{1};
    end
    
    if nargin <= 8
        K = eye(3);
        W = 1920;W_1 = W - 1;
        H = 1080;H_1 = H - 1;
    else
        K = varargin{2};
        W = round(K(1,3)*2);W_1 = W - 1;
        H = round(K(2,3)*2);H_1 = H - 1;
    end
 
    if nargin <= 9
        isplanar = 0;
    else
        isplanar = varargin{3};
    end
    
    d1 = abs(dot(nw,C1)+dw);  
    miu = refra_n2/refra_n1;%
    miu = 1/miu;
    for i = 1:N        
        while true
        f1 = [rand(1)*W_1+1;rand(1)*H_1+1;1];x1 = inv(K)*f1;
        x1 = x1 ./ norm(x1);
        q1 = R1'*x1;q1 = q1 ./ norm(q1);
        ref_A = -d1 / dot(q1,nw) * q1 + C1;
        i1 = q1; % incident ray 1
        nn = -nw;
        tr1 = miu * i1 + sqrt(1-miu^2*(1-(dot(nn,i1))^2)) * nn - miu * dot(nn,i1) * nn;% transmitted ray 1
        tr1 = tr1./norm(tr1);
        if ~isreal(tr1)
            continue;
        end
        
        if ~isplanar
            Q = ref_A + tr1 * (( 0.5 + 0.5 * rand(1) )/abs(tr1(3)));%+ 0.3 * rand(1) %+ 0.5 * rand(1) + 0.5 * rand(1)
        else
            aa = -ref_A(3) / tr1(3);
            % planar
            Q = ref_A + tr1 * aa;
        end
        
        break
        end
        f1(1:2) = f1(1:2) + std.*randn(2,1);%q1 = q1 ./ norm(q1);
        q1s(:,i) = f1;
        Q1s(:,i) = Q;
        ref_As(:,i) = ref_A;
        plane_norm = cross(tr1,nn);plane_norm = plane_norm ./ norm(plane_norm);
        d = -dot(ref_A,plane_norm);
    end
    varargout{1} = q1s(1:2,:);
%     figure
%     plot(q1s(1,:),q1s(2,:),'r.');
%     xlim([1,1920]);
%     ylim([1,1080]);
    varargout{2} = Q1s;
    varargout{3} = ref_As;
end