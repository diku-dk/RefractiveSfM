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

function [Cost] = CostFuncOneLayersCase1(R,C,n,d,rayDnoise,XYZ,KK,mu)
    N = size(XYZ,2);
    Cost = zeros(2*N,1);
    ps = R'*(XYZ - C);
    for ii = 1:N
        p = ps(:,ii);
        M = SolveFlatRefractionPointOneLayerCase1(d,n,mu,p);
        M = KK*M;
        px1 = M(1)/M(3);
        py1 = M(2)/M(3);
        px = rayDnoise(1,ii);
        py = rayDnoise(2,ii);
        Cost(2*ii-1,1) = px-px1;
        Cost(2*ii,1) = py-py1;
    end
end