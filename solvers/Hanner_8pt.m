function [R,t] = Hanner_8pt(Q,q,K,miu,varargin)
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
    R1 = varargin{1};

    %% use cvpr15 solution
    if size(q,1) == 2
        q = [q;ones(1,size(q,2))];
    end
    q = inv(K)*q;
    v = pose3dsolver_5pt(Q(:,:),q(1:2,:),'r',1/miu);

    minerr = 1e6;
    minid = -1;
    for jj = 1:size(v,2)
        Rhat = (q2R(v(1:4,jj)));
        err = norm(R1'*Rhat-eye(3),'fro');
        if err < minerr
            minerr = err;
            minid = jj;
        end
    end
    if isempty(v)
        R = eye(3);
        t = zeros(3,1);
    else
        R = q2R(v(1:4,minid));
        t = v(5:end,minid);
    end  
end