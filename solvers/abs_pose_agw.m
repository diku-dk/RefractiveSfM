function [R,C] = abs_pose_agw(Q,q,K,ior1,ior2,varargin)
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
% varargin:
%   1: run iterative optimization or not.
%   2: is the planar case.
%   3: verbose.
    apply_iter_opt = 0;
    verbose = 1;
    is_plane = 0;
    if ~isempty(varargin)
        apply_iter_opt = varargin{1};
        if length(varargin) > 1
            is_plane = varargin{2};
        end
        if length(varargin) > 2
            verbose = varargin{3};
        end
    end
    mu = ior2/ior1 ; % refractive index of medium

    if size(q,1) == 2
        q = [q;ones(1,size(q,2))];
    end
    Pall = Q;
    ray2D = inv(K)*q;
    
    if verbose == 1
        disp('=================================')
        disp('=================================')
        disp('==========8 pt =======================')
    end
    
    if is_plane == 0
        solVec = EightPointEstimationOneLayerCase1KnownMu(ray2D,Pall,mu,300,0);
    else
        shift = [0;0;-mean(Pall(3,:))];
        Pall = Pall+shift;
        solVec = AxisEstimationPlanarSceneOneLayerCase1(ray2D,Pall,mu);
    end
    ne8pt = solVec(1:3);
    tau8pt = solVec(4);
    R8pt = RPY2Rot(solVec(5:7));
    t8pt = solVec(8:10);
    w8pt = solVec(5:7);

    x0 = [ne8pt(1:2)/abs(ne8pt(3)) ; w8pt ; t8pt ; tau8pt ; ];
    if apply_iter_opt == 1
        if verbose == 1
            disp('===============================')
            disp('===============================')
            disp('===============================')
            fprintf('We got the initial estimates. Now minimizing re-projection error using analytical forward projection equations');
        end
        %non-linear optimization
        options = optimset('Display','off','TolFun',1e-6,'TolX',1e-6,'MaxIter',5000,'MaxFunEvals',10000);
        x = lsqnonlin(@CostFuncOneLayersCase1,x0,[],[],options,q,Pall,K,mu);
    else
        x = x0;
    end
    R = RPY2Rot(x(3:5));
    t = x(6:8);
    C = -R'*t;
    if is_plane == 1
        C = C - shift;
    end
end
