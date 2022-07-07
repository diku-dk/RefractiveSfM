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



function [vr,a,b] = RefractedRay(vi,n,mu1,mu2)
    % vi is incoming ray
    % n is normal
    % mu1 and mu2 are refractive indices at the boundary
    % computes the outgoing refracted ray
        
    n = n/norm(n);

    a = mu1/mu2;
    b = -mu1*(vi'*n) - sqrt(mu1^2*(vi'*n)^2 - (mu1^2-mu2^2)*(vi'*vi));
    b = b/mu2;

    vr = a*vi + b*n;





