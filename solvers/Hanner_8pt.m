function [R,t] = Hanner_8pt(Q,q,K,miu,varargin)
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