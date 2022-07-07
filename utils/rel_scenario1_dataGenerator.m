%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright 2022, Xiao Hu, Francois Lauze, Kim Steenstrup Pedersen
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License. 
%   You may obtain a copy of the License at
%
%       http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = dataGenerator(N,R1,R2,n1,n2,t1,t2,refra_n1,refra_n2,ds1,ds2,varargin)
% R1, R2 from: w to c
    if nargin < 12
        std = 0*1e-6;
    else
        std = varargin{1};
    end
    if nargin < 13
        KK = eye(3);
    else
        KK = varargin{2};
    end
    
    d1 = abs(dot(n1,t1)+ds1);  d2 = abs(dot(n2,t2)+ds2);
    miu = refra_n2/refra_n1;%
    
    w=KK(1,3)*2;
    h=KK(2,3)*2;
    
    for i = 1:N        
        while true
        f1 = [rand(1)*(w-1)+1;rand(1)*(h-1)+1;1];x1 = inv(KK)*f1;% ray in c
%         x1 = [rand(2,1)*4-2;1];x1 = x1./norm(x1);
        q1 = R1'*x1;q1 = q1 ./ norm(q1);% ray in w
        ref_A = d1 / dot(q1,n1) * q1 + t1;% refractive point in w
        i1 = q1; % incident ray 1
        nn = n1;
        tr1 = miu * i1 + sqrt(1-miu^2*(1-(dot(nn,i1))^2)) * nn - miu * dot(nn,i1) * nn;% transmitted ray 1

        if ~isreal(tr1)
            continue;
        end
        
        tr1 = tr1./norm(tr1);
        Q = ref_A + tr1 * (rand(1)*5+5);% world point
        [q2,~,ref_B] = forward_project(n2,d2,refra_n1,refra_n2,t2,Q);% refractive point in another plane
    
        i2 = q2; % incident ray 1
        nn = n2;
        tr2 = miu * i2 + sqrt(1-miu^2*(1-(dot(nn,i2))^2)) * nn - miu * dot(nn,i2) * nn;% incident ray 2
        
        if ~isreal(i2)
            continue;
        end
        if ~isreal(tr2)
            continue;
        end
        
        l1 = plucker_line([ref_A;1],[tr1;0]);
        l2 = plucker_line([ref_B;1],[tr2;0]);

        if norm(cross(n2,[0,0,1]')) < 1e-6
            assert(norm(-t2(3)/q2(3)*q2+t2 - ref_B)<1e-6)
        else
            assert(norm(t2+d2/dot(q2,n2)*q2-ref_B)<1e-6)
        end
        assert(norm(cross(tr1,Q-ref_A))<1e-6)
        assert(norm(cross(tr2,Q-ref_B))<1e-6)
        assert(norm(l1([4,5,6,1,2,3])'*l2)<1e-6)
        assert(norm(l2(4:6)+cross(l2(1:3),Q))<1e-6)
        break
        end
        % to camera 2
        q1 = R1 * q1;% ray in camera
        q2 = R2 * q2;% ray in camera
%         q1 = q1 + std.*randn(3,1);q1 = q1 ./ norm(q1);
%         q2 = q2 + std.*randn(3,1);q2 = q2 ./ norm(q2);
        imgf1 = KK * q1;imgf1 = imgf1 ./ imgf1(3,:);
        imgf2 = KK * q2;imgf2 = imgf2 ./ imgf2(3,:);
        q1s(:,i) = imgf1 + [std.*randn(2,1);0];
        q2s(:,i) = imgf2 + [std.*randn(2,1);0];
        ref_As(:,i) = ref_A;
        ref_Bs(:,i) = ref_B;
        Qs(:,i) = Q;
    end
    varargout{1} = q1s;
    varargout{2} = q2s;
    varargout{3} = ref_As;
    varargout{4} = ref_Bs;
    varargout{5} = Qs;
end