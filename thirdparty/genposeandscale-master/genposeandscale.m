function [Psolns,ssolns] = genposeandscale(p,d,q)
% GENPOSEANDSCALE Solves for the pose and scale of a generalized camera.
%   [Psolns,ssolns] = genposeandscale(p,d,q)
%
%   p is a 3xN matrix of ray origins
%   d is a 3xN matrix of ray directions
%   q is a 3xN matrix of anchor points
%
%   N must be at least 4.
%
% If you use this algorithm in your work, please cite the following paper:
% J. Ventura, C. Arth, G. Reitmayr, and D. Schmalstieg,
% "A Minimal Solution to the Generalized Pose-and-Scale Problem", CVPR, 2014.
% 
% Copyright (c) 2014 Graz University of Technology.
% 
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted (subject to the limitations in the disclaimer
% below) provided that the following conditions are met:
% 
%  * Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 
%  * Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
% 
%  * Neither the name of Graz University of Technology nor the names of its
%    contributors may be used to endorse or promote products derived from this
%    software without specific prior written permission.
% 
% NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
% LICENSE.  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

B = getbasis(p,d,q);
[b1,b2,b3,b4,b5] = genposeandscale_solvecoeffs(B(:,1),B(:,2),B(:,3),B(:,4),B(:,5),B(:,6));

nsolns = length(b1);
Psolns = cell(nsolns,1);
ssolns = zeros(nsolns,1);
for i=1:nsolns,
    b = [b1(i);b2(i);b3(i);b4(i);b5(i);1];
    P = B*b;
    
    myr=P(1:9);
    myR=reshape(myr,3,3);
    
    a = norm(myR(1,:));
    
    if det(myR)<0,
        a=-a;
    end
    
    myt=P(10:12);
    mys=P(13);
    
    myR = myR/a;
    myt = myt/a;
    mys = mys/a;
    
    Psolns{i}=[myR' -myR'*myt];
    ssolns(i)=1/mys;
end

function B = getbasis(p,d,q)
N = size(p,2);
A = zeros(N*2,13);
n = 1;
for i=1:N,
	A(n:n+1,:) = getcross(p(:,i),d(:,i),q(:,i));
	n = n+2;
end
[~,~,V] = svd(A);
B = V(:,8:13);

function A = getcross(p,d,q)
p11 = p(1);
p12 = p(2);
p13 = p(3);

d11 = d(1);
d12 = d(2);
d13 = d(3);

q11 = q(1);
q12 = q(2);
q13 = q(3);

A = [
	d13*q11,       0, -d11*q11, d13*q12,       0, -d11*q12, d13*q13,       0, -d11*q13, d13,   0, -d11, d11*p13 - d13*p11;
          0, d13*q11, -d12*q11,       0, d13*q12, -d12*q12,       0, d13*q13, -d12*q13,   0, d13, -d12, d12*p13 - d13*p12
];
