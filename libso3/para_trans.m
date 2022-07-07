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

function vpt = para_trans(R1,R2,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     A = logSO3(R1'*R2);
%     Rpt = expSO3(A/2);
%     vpt = invhat(R2'*R1*Rpt*hat(v)*Rpt);
    % fast version
    vpt = expSO3(logSO3(R2'*R1)/2)*v;
end


%function vpt = para_trans(R1,R2,v)
%    A = logSO3(R1'*R2);
%    Rpt = expSO3(A/2);
%    if size(v,2) == 1
%        vpt = invhat(R2'*R1*Rpt*hat(v)*Rpt);
%    else
%        vpt = invhat(R2'*R1*Rpt*(R1'*v)*Rpt);
%    end
%end