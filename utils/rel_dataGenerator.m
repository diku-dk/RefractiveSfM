function varargout = rel_dataGenerator(N,R1,n1w,d1w,C1,C2,refra_n1,refra_n2,varargin)
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
%   R1: True camera rotation from the camera 1 frame to world frame.
%   C1,C2: Camera 1 and 2's location in thr world frame.
%   n1w,d1w: parameters of refractive interface.
%   refra_n1,refra_n2: refractive indices.
% Optional:
%   std: Noise standard deviation.
%   K: Camera Intrinsic Matrix.
%
% Outputs:
%   q1s,q2s: feature points.
%   Q1s: 3D points.
%   ref_As,ref_Bs: refractive points.
%   R2: True camera rotation from the camera 2 frame to world frame.
    
    % R1, R2 from: c to w
    if nargin < 9
        std = 0*1e-6;
    else
        std = varargin{1};
    end
    if nargin < 10
        KK = eye(3);
    else
        KK = varargin{2};
    end
    
    d1 = abs(dot(n1w,C1)+d1w);  d2 = abs(dot(n1w,C2)+d1w);
    miu = refra_n2/refra_n1;%
    
    w=KK(1,3)*2;
    h=KK(2,3)*2;
    
    q2s = zeros(3,N);
    for i = 1:N        
        while true
        f1 = [rand(1)*(w-1)+1;rand(1)*(h-1)+1;1];x1 = inv(KK)*f1;% ray in c
        q1 = R1*x1;q1 = q1 ./ norm(q1);% ray in w
        ref_A = d1 / dot(q1,n1w) * q1 + C1;% refractive point in w
        i1 = q1; % incident ray 1
        nn = n1w;
        tr1 = miu * i1 + sqrt(1-miu^2*(1-(dot(nn,i1))^2)) * nn - miu * dot(nn,i1) * nn;% transmitted ray 1

        if ~isreal(tr1)
            continue;
        end
        
        tr1 = tr1./norm(tr1);
        Q = ref_A + tr1 * (rand(1)*5+5);% world point
        [q2,~,ref_B] = forward_project(n1w,d2,refra_n1,refra_n2,C2,Q);% refractive point in another plane
    
        i2 = q2; % incident ray 1
        nn = n1w;
        tr2 = miu * i2 + sqrt(1-miu^2*(1-(dot(nn,i2))^2)) * nn - miu * dot(nn,i2) * nn;% incident ray 2
        
        if ~isreal(i2)
            continue;
        end
        if ~isreal(tr2)
            continue;
        end
        
        l1 = plucker_line([ref_A;1],[tr1;0]);
        l2 = plucker_line([ref_B;1],[tr2;0]);

        if norm(cross(n1w,[0,0,1]')) < 1e-6 && abs(d2) < 1e-6
            assert(norm(-C2(3)/q2(3)*q2+C2 - ref_B)<1e-6)
        else
            assert(norm(C2+d2/dot(q2,n1w)*q2-ref_B)<1e-6)
        end
        assert(norm(cross(tr1,Q-ref_A))<1e-6)
        assert(norm(cross(tr2,Q-ref_B))<1e-6)
        assert(norm(l1([4,5,6,1,2,3])'*l2)<1e-6)
        assert(norm(l2(4:6)+cross(l2(1:3),Q))<1e-6)
        break
        end
        % to camera 2
        q1 = R1' * q1;% ray in camera
        q2s(:,i) = q2;
        
        imgf1 = KK * q1;imgf1 = imgf1 ./ imgf1(3,:);
        q1s(:,i) = imgf1 + [std.*randn(2,1);0];
        ref_As(:,i) = ref_A;
        ref_Bs(:,i) = ref_B;
        Qs(:,i) = Q;
    end
    
    % find a R2 that best covers all q2
    mean_ref_B = mean(ref_Bs,2);
    target_ray = (mean_ref_B - C2);
    target_ray = target_ray ./ norm(target_ray);
    rotaxis = cross([0;0;1],target_ray);
    rotaxis = rotaxis ./ norm(rotaxis);
    theta = acos(dot([0;0;1],target_ray));
    R = expSO3(theta*rotaxis);
    assert(abs(dot(R*[0;0;1],target_ray)-1)<1e-6);% from c to w
    R2 = R;% from c to w
    
    for i = 1:size(q2s,2)
        q2 = q2s(:,i);
        q2 = R2' * q2;% ray in camera
        imgf2 = KK * q2;imgf2 = imgf2 ./ imgf2(3,:);
        q2s(:,i) = imgf2 + [std.*randn(2,1);0];
    end
    
    varargout{1} = q1s;
    varargout{2} = q2s;
    varargout{3} = ref_As;
    varargout{4} = ref_Bs;
    varargout{5} = Qs;
    varargout{6} = R2;
end
