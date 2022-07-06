function [R,C] = abs_pose_agw(Q,q,K,ior1,ior2,varargin)
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
