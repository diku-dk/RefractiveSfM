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

function [X,R,C,timelsq] = refinement_chadebecq_iccv(imgf1s,imgf2s,n1,n2,d1,d2,ior1,ior2,miu,KK,R,C)
    q1s = imagepoints_to_rays(KK,imgf1s);
    q2s = imagepoints_to_rays(KK,imgf2s);

    r11 = imagepoints_to_rays(eye(3),double(q1s));
    r12 = compute_refractive_ray(ior1,ior2,r11,n1);
    p1 = compute_refractive_point(r11,d1,n1);
    r21 = imagepoints_to_rays(eye(3),double(q2s));
    r22 = compute_refractive_ray(ior1,ior2,r21,n2);
    p2 = compute_refractive_point(r21,d2,n2);

    lparams1 = zeros(6,size(p1,2));
    for i = 1:size(p1,2)
        lparams1(:,i) = plucker_line([p1(:,i);1],[r12(:,i);0]);
    end
    lparams2 = zeros(6,size(p2,2));
    for i = 1:size(p2,2)
        lparams2(:,i) = plucker_line([p2(:,i);1],[r22(:,i);0]);
    end
    
    % triangulation under refraction
    ray1 = r12; Q1 = p1;
    ray2 = R * r22; Q2 = R * p2 + C;
    X = zeros(3,size(imgf1s,2));
    for ii = 1:size(Q2,2)
        AA = [eye(3) -ray1(1:3,ii) zeros(3,1);eye(3) zeros(3,1) -ray2(1:3,ii)];
        bb = [Q1(1:3,ii);Q2(1:3,ii)];
        sol = AA\bb;
        X(:,ii) = sol(1:3);                
    end   

    %non-linear optimization
    x0 = [X(:); logSO3(R) ; C];x0 = double(x0);
    options = optimset('Display','off','TolFun',1e-6,'TolX',1e-6,'MaxIter',1000,'MaxFunEvals',10000,'Algorithm','levenberg-marquardt');
    tic;
    x = lsqnonlin(@cost_func_total,x0,[],[],options,imgf1s,imgf2s,n1,n2,d1,d2,lparams1,lparams2,miu,KK);
    timelsq = toc;
    X = reshape(x(1:end-6),3,[]);    
    R = expSO3(x(end-5:end-3));
    C = x(end-2:end);
end

function [cost] = cost_epipolar(R,C,l1,l2)
    tmp1 = R * l2(1:3,:);
    tmp2 = R * l2(4:6,:);
    Cx = skewm(C);
    err = dot(l1(1:3,:),Cx*tmp1) + dot(l1(4:6,:),tmp1) + dot(l1(1:3,:),tmp2);
    cost = err';
end

function cost = cost_func_total(x,q1,q2,n1,n2,d1,d2,l1,l2,miu,KK)
    X = reshape(x(1:end-6),3,[]);    
    R = expSO3(x(end-5:end-3));
    C = x(end-2:end);

    cost1 = CostFuncOneLayersCase1(eye(3),zeros(3,1),-n1,d1,q1,X,KK,miu);
    cost2 = CostFuncOneLayersCase1(R,C,-n2,d2,q2,X,KK,miu);
    cost3 = cost_epipolar(R,C,l1,l2);
    
    cost = [cost1;cost2;cost3];
end

