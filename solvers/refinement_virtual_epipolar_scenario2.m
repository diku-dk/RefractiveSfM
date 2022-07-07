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


function [X,R,C,timelsq] = refinement_virtual_epipolar_scenario2(n1w, d1w, ior1, ior2, imgf1s, imgf2s, KK, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code for relative pose estimation algorithm proposed in:
%
%       Absolute and Relative Pose Estimation in Refractive Multi View
%       In Submission to ICCV 2021
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% new section: relative pose refineemnt.
    rotaxis = cross(n1w,[0;0;1]);
    rotaxis = rotaxis ./ norm(rotaxis);
    theta = acos(dot(n1w,[0;0;1]));
    if abs(theta) < 1e-6
        Rpre = eye(3);
    else
        Rpre = expSO3(theta*rotaxis);
        assert(abs(dot(Rpre*n1w,[0;0;1])-1)<1e-6);
    end
    Cpre = -d1w*n1w;
    
%     Qpre = Rpre * (Qs-Cpre);
%     n1pre = [0;0;1];
%     d1pre = 0;
    
    lambda = ior1/ior2;


    % plot to check
%     C1pre = Rpre * (COP1-Cpre);
%     C2pre = Rpre * (COP2-Cpre);
% 
%     ref_As_pre = Rpre * (ref_As-Cpre);
%     ref_Bs_pre = Rpre * (ref_Bs-Cpre);
%     figure;cmap = lines(9);
%     plot3(0,0,0,'m-s','MarkerSize',15);hold on;
%     plot3(C1pre(1),C1pre(2),C1pre(3),'-d','Color',cmap(1,:),'MarkerSize',15);hold on;
%     h1=plotCameraFrustum(R1'*Rpre',-(R1'*Rpre')*C1pre,KK,0.5,cmap(1,:),1.5);
%     plot3(C2pre(1),C2pre(2),C2pre(3),'-d','Color',cmap(2,:),'MarkerSize',15);hold on;
%     plotCameraFrustum(R2'*Rpre',-(R2'*Rpre')*C2pre,KK,0.5,cmap(2,:),1.5);
% 
%     [xx,yy] = meshgrid(min(ref_As_pre(1,:))-1:0.1:max(ref_As_pre(1,:))+1,min(ref_As_pre(2,:))-1:0.1:max(ref_As_pre(2,:))+1);
%     zz = (-d1pre - n1pre(1)*xx - n1pre(2)*yy)/n1pre(3);
%     h3=mesh(xx,yy,zz);hold on
%     for ii = 1:num_features
%         h=line([C1pre(1) ; ref_As_pre(1,ii)], [C1pre(2) ; ref_As_pre(2,ii)],[C1pre(3) ; ref_As_pre(3,ii)] ,'LineWidth', 2 , 'Color', cmap(3,:));hold on;
%         h.Color(4) = 0.5;
%         h2=plot3(ref_As_pre(1,ii),ref_As_pre(2,ii),ref_As_pre(3,ii),'s', 'Color',cmap(4,:),'MarkerFaceColor','y');hold on;
%         h=line([Qpre(1,ii) ; ref_As_pre(1,ii)], [Qpre(2,ii) ; ref_As_pre(2,ii)],[Qpre(3,ii) ; ref_As_pre(3,ii)] ,'LineWidth', 2, 'Color', cmap(5,:));
%         h.Color(4) = 0.5;
% 
%         h=line([C2pre(1) ; ref_Bs_pre(1,ii)], [C2pre(2) ; ref_Bs_pre(2,ii)],[C2pre(3) ; ref_Bs_pre(3,ii)] ,'LineWidth', 2 , 'Color', cmap(6,:));hold on;
%         h.Color(4) = 0.5;
%         plot3(ref_Bs_pre(1,ii),ref_Bs_pre(2,ii),ref_Bs_pre(3,ii),'o', 'Color',cmap(7,:),'MarkerFaceColor','y');hold on;
%         h=line([Qpre(1,ii) ; ref_Bs_pre(1,ii)], [Qpre(2,ii) ; ref_Bs_pre(2,ii)],[Qpre(3,ii) ; ref_Bs_pre(3,ii)] ,'LineWidth', 2, 'Color', cmap(8,:));
%         h.Color(4) = 0.5;
%     end
%     h4=plot3(Qpre(1,:),Qpre(2,:),Qpre(3,:),'ko','MarkerFaceColor','k');hold on;
%     xlabel('x: (m)','Interpreter','latex');ylabel('y: (m)','Interpreter','latex');zlabel('z: (m)','Interpreter','latex');
%     title('Refraction in idea world coordinate','Interpreter','latex');
%     legend([h1,h2,h3,h4],{'Camera idea','Refractive Point idea','Refractive Plane idea','World Point idea'},'Interpreter','latex');
    
    %% solution-1: use normal essential matrix
    if isempty(varargin)
        E = cv.findEssentialMat(imgf1s(1:2,:)', imgf2s(1:2,:)', 'CameraMatrix', KK);
        [R, t, good, ~] = cv.recoverPose(E, imgf1s(1:2,:)', imgf2s(1:2,:)', ...
            'CameraMatrix', KK, 'DistanceThreshold', 50);
        C0 = -R'*t;
        R0 = R';
    else
        R0 = varargin{1};
        C0 = varargin{2};
    end
    [X,R,C,timelsq] = refinement_vc(imgf1s,imgf2s,ior1,ior2,lambda,KK,Rpre,Cpre,R0,C0);
    
    
end

function [X,R,C,timelsq] = refinement_vc(f1s,f2s,ior1,ior2,miu,K,Rpre,Cpre,R0,C0)
    r11 = imagepoints_to_rays(K,f1s);% K-1 u  and normalize
    r21 = imagepoints_to_rays(K,f2s);% K-1 u  and normalize
        
    % C in idea world
    C11 = -Rpre * Cpre;
    C21 = Rpre * (C0 - Cpre);
    
    if C11(3) < 0
        n_in_iw = [0 0 1]';
    else
        n_in_iw = [0 0 -1]';
    end
    
    % to world
    r12 = r11;
    r22 = R0 * r21;
    
    % to idea world
    r13 = Rpre * r12;
    r23 = Rpre * r22;
    
    % check
    r14 = [miu.*r13(1,:);miu.*r13(2,:);n_in_iw(3).*sqrt(1-miu^2+miu^2*r13(3,:).^2)];
    r24 = [miu.*r23(1,:);miu.*r23(2,:);n_in_iw(3).*sqrt(1-miu^2+miu^2*r23(3,:).^2)];
    
    % compute refractive point
    q11 = compute_refractive_point(r13,n_in_iw,C11);
    q21 = compute_refractive_point(r23,n_in_iw,C21);
    
    % triangulation under refraction
    ray1 = r14; Q1 = q11;
    ray2 = r24; Q2 = q21;
    X = zeros(3,size(f1s,2));
    for ii = 1:size(Q2,2)
        AA = [eye(3) -ray1(1:3,ii) zeros(3,1);eye(3) zeros(3,1) -ray2(1:3,ii)];
        bb = [Q1(1:3,ii);Q2(1:3,ii)];
        sol = AA\bb;
        X(:,ii) = sol(1:3);                
    end   
    X = Rpre' * X + Cpre;
    
    x0 = [X(:); logSO3(R0) ; C0];
    x0 = double(x0);
    tic;
    x = solver_levmar_own(x0,miu,K,Rpre,Cpre,r11,r21);
    timelsq = toc;

    %non-linear optimization
    X = reshape(x(1:end-6),3,[]);    
    R = expSO3(x(end-5:end-3));
    C = x(end-2:end);
end

function x = solver_levmar_own(x,lambda,K,Rpre,Cpre,r11,r21)
    k = 0;
    v = 2;
    
    J = [];
    err = [];
    
    [err,~,~,J] = cost_func_total_own(x,lambda,K,Rpre,Cpre,r11,r21);
    A = sparse(J'*J); g = sparse(J'*err);
    errTol = 1e-10;
    found = max(abs(g)) < errTol;
    tau = 1e-3;
    miu = tau * max(diag(A));
    maxIter = 1000;
    iter = 0;
    epsilon2 = 1e-10;
     % routine for lm
    disp(['init cost is:',num2str((err'*err))]);
    while found == false && iter < maxIter
        iter = iter + 1;
        hlm = (A+miu*speye(size(A,1),size(A,2)))\(-g);
        hnorm = sqrt(hlm'*hlm);
        if hnorm <= epsilon2*(sqrt(x'*x)+epsilon2)
            found = true;
        else
            X = reshape(x(1:end-6),3,[]);
            R = expSO3(x(end-5:end-3));
            C = x(end-2:end);
            newX = X + reshape(hlm(1:end-6),3,[]);
            newR = expSO3(logSO3(R))*expSO3(hlm(end-5:end-3));
            newC = C + hlm(end-2:end);
            newx = [newX(:); logSO3(newR) ; newC];
            [newerr,~,~,newJ] = cost_func_total_own(newx,lambda,K,Rpre,Cpre,r11,r21);
            err_gain = (err'*err) - (newerr'*newerr);
            gain = err_gain / (0.5*(hlm'*(miu*hlm-g)));
            if gain > 0
                disp(['iter:',num2str(iter),' cost is:',num2str((newerr'*newerr))]);
                x = newx;
                err = newerr;
                A = sparse(newJ'*newJ);
                g = sparse(newJ'*newerr);
                found = max(abs(err)) < errTol;
                miu = miu * max(1/3,1-(2*gain-1)^3);v = 2;
                oldcost = sqrt(newerr'*newerr);
            else
                miu = miu * v;v = v*2;
            end
        end
    end
    disp(['final cost is:',num2str((newerr'*newerr))]);
end

function varargout = cost_func_total_own(x,miu,K,Rpre,Cpre,r11,r21)
    X0 = reshape(x(1:end-6),3,[]);
    R0 = expSO3(x(end-5:end-3));
    C0 = x(end-2:end);

    % C in idea world
    C11 = -Rpre * Cpre;
    C21 = Rpre * (C0 - Cpre);
    
    if C11(3) < 0
        n_in_iw = [0 0 1]';
    else
        n_in_iw = [0 0 -1]';
    end
    
    % to world
    r12 = r11;
    r22 = R0 * r21;
    
    % to idea world
    r13 = Rpre * r12;
    r23 = Rpre * r22;
    
    % check
    r14 = [miu.*r13(1,:);miu.*r13(2,:);n_in_iw(3).*sqrt(1-miu^2+miu^2*r13(3,:).^2)];
    r24 = [miu.*r23(1,:);miu.*r23(2,:);n_in_iw(3).*sqrt(1-miu^2+miu^2*r23(3,:).^2)];
    
    % compute refractive point
    q11 = [C11(1)*r13(3,:)-C11(3)*r13(1,:);...
            C11(2)*r13(3,:)-C11(3)*r13(2,:);...
            zeros(1,size(r13,2))];
    q11 = q11 ./ r13(3,:);
    
    q21 = [C21(1)*r23(3,:)-C21(3)*r23(1,:);...
            C21(2)*r23(3,:)-C21(3)*r23(2,:);...
            zeros(1,size(r23,2))];
    q21 = q21 ./ r23(3,:);
    
    % compute virtual camera transformation
    [~,Cv1] = find_virtual_camera_coordiante(q11,r14,n_in_iw,C11,miu,r13);
    [~,Cv2] = find_virtual_camera_coordiante(q21,r24,n_in_iw,C21,miu,r23);
    
    % to virtual camera 
    q12 = to_virtual_camera_coordinate(q11,[],[],eye(3),Cv1);
    q22 = to_virtual_camera_coordinate(q21,[],[],eye(3),Cv2);
    
    [err1,AtA1,Atb1,Js1] = cost_func_reproj(X0,eye(3),Cv1,K,q12,Rpre,Cpre,1);
    [err2,AtA2,Atb2,Js2] = cost_func_reproj(X0,eye(3),Cv2,K,q22,Rpre,Cpre,0,n_in_iw,miu,r21,r23,r24,C21,R0);
    [err3,AtA3,Atb3,Js3] = cost_epipolar_own(r14,r24,Cv1,Cv2,miu,n_in_iw,Rpre,R0,r21,r23,C21);
    
    AtA = AtA1 + AtA2 + AtA3;%
    Atb = Atb1 + Atb2 + Atb3;%
    err = [err1;err2;err3];%err1;
    Js = [Js1;Js2;Js3];
    varargout{1}=err;
    varargout{2}=AtA;
    varargout{3}=Atb;
    varargout{4}=Js;
end


function varargout = cost_func_reproj(X,Rv1,Cv1,K,qv,Rpre,Cpre,omitRt,varargin)
    if omitRt == 0
        n_in_iw = varargin{1};
        miu = varargin{2};
        r1 = varargin{3};
        r3 = varargin{4};
        r4 = varargin{5};
        C1 = varargin{6};
        R0 = varargin{7};
    end

    N = size(X,2);
    M = 3*N + 6;
    fv = K(1,1);
    Cvs = horzcat(Cv1{:});
    Xv = Rv1' * (Rpre * (X - Cpre)) - Rv1' * Cvs;
    gv = [fv.*Xv(1,:)./Xv(3,:) - fv.*qv(1,:)./qv(3,:);...
          fv.*Xv(2,:)./Xv(3,:) - fv.*qv(2,:)./qv(3,:)];
    err = gv(:);
    
    fvXv3 = fv./Xv(3,:);
    fvXv1Xv3 = -fv.*Xv(1,:)./(Xv(3,:).^2);
    fvXv2Xv3 = -fv.*Xv(2,:)./(Xv(3,:).^2);
    
    fvq2v3 = -fv./qv(3,:);
    fvq21Xv3 = fv.*qv(1,:)./(qv(3,:).^2);
    fvq22Xv3 = fv.*qv(2,:)./(qv(3,:).^2);
    
    AtA = zeros(M,M);
    Atb = zeros(M,1);
    
    Js = zeros(length(err),M);
    
    for i = 1:N
        J = zeros(2,M);
        
        jacXv = [fvXv3(i) 0 fvXv1Xv3(i);...
                0 fvXv3(i) fvXv2Xv3(i)];
        jacqv = [fvq2v3(i) 0 fvq21Xv3(i);...
                0 fvq2v3(i) fvq22Xv3(i)];
            
        J(:,i*3-2:i*3) = jacXv * Rpre;
        
        if omitRt == 0
            J3 = -eye(3);
            J4 = [0 0 0;0 0 0;0 0 -C1(3)*n_in_iw(3)*(1-miu^2)/(miu*r3(3,i)^2*r4(3,i)/n_in_iw(3))] * Rpre * (-R0 * skewm(r1(:,i)));
            JvvR = J3 * J4;
        
            J6 = [1 0 0;0 1 0;0 0 r4(3,i)/miu/r3(3,i)] * Rpre;
            JvvC = J3 * J6;
            
            J7 = eye(3);
            JqR = J7 * [-C1(3)/r3(3,i) 0 C1(3)*r3(1,i)/(r3(3,i)^2); ...
                   0 -C1(3)/r3(3,i) C1(3)*r3(2,i)/(r3(3,i)^2); 0 0 0]  * Rpre * (-R0 * skewm(r1(:,i)));
                
            JqC = J7 * [1 0 -r3(1,i)/r3(3,i); ...
                   0 1 -r3(2,i)/r3(3,i); 0 0 0]  * Rpre;
            
            J(:,end-5:end-3) = jacXv * JvvR + jacqv * JqR + jacqv * JvvR;
            J(:,end-2:end) = jacXv * JvvC + jacqv * JqC + jacqv * JvvC;
        end
        AtA = AtA + J'*J;
        Atb = Atb + J'*err(i*2-1:i*2);
        
        Js(i*2-1:i*2,:) = J;
    end
    varargout{1}=err;
    varargout{2}=AtA;
    varargout{3}=Atb;
    varargout{4}=Js;
end

function varargout = cost_epipolar_own(r14,r24,Cv1,Cv2,miu,n_in_iw,Rpre,R0,r21,r23,C21)
    N = size(r14,2);
    M = 3*N + 6;    
    err = zeros(N,1);
    AtA = zeros(M,M);
    Atb = zeros(M,1);
    
    Js = zeros(length(err),M);
    for i = 1:size(r14,2)
        J = zeros(1,M);
        rF = skewm(Cv2{i}-Cv1{i});
        err(i) = r24(:,i)'*rF*r14(:,i);
        
        J1 = (rF*r14(:,i))';
        J2 = [miu,0,0;0,miu,0;0 0 n_in_iw(3)*n_in_iw(3)/r24(3,i)*miu^2*r23(3,i)] * Rpre * (-R0 * skewm(r21(:,i)));
        Jr24R = J1 * J2;
        
        J3 = -r24(:,i)'*skewm(r14(:,i));
        J4 = [0 0 0;0 0 0;0 0 -C21(3)*n_in_iw(3)*(1-miu^2)/(miu*r23(3,i)^2*r24(3,i)/n_in_iw(3))] * Rpre * (-R0 * skewm(r21(:,i)));
        JvvR = J3 * J4;
        
        J6 = [1 0 0;0 1 0;0 0 r24(3,i)/miu/r23(3,i)] * Rpre;
        JvvC = J3 * J6;
        
        J(:,end-5:end) = [Jr24R+JvvR JvvC];
        AtA = AtA + J'*J;
        Atb = Atb + J'*err(i);
        Js(i,:) = J;
    end
    varargout{1} = err;
    varargout{2}=AtA;
    varargout{3}=Atb;
    varargout{4}=Js;
end

function [Rv,Cv,varargout] = find_virtual_camera_coordiante(ps,rws,n,C,lambda,r3)
% this function computes the virtual camera coordiante system for RSFM as
% described in the paper:
% Refractive Structure-from-Motion on Underwater Images
    n = n./norm(n);
    rotaxis = cross([0;0;1],n);
    rotaxis = rotaxis./norm(rotaxis);
    rotangle = acos(dot([0;0;1],n));
    if abs(rotangle) < 1e-6
        Rv = eye(3);
    else
        Rv = expSO3(rotaxis*rotangle)';
    end
        
    Cv = cell(size(rws,2),1);
    tv = cell(size(rws,2),1);
    for i = 1:size(rws,2)
        rw = rws(:,i); p = ps(:,i);
        rw = rw./norm(rw);
        A = [n -rw];b = p - C;
        res = A\b;
        
        lambda1 = C(3)/r3(3,i)/lambda;
        assert(abs(res(2)-lambda1)<1e-8);
        
        Cv{i} = p + res(2)*rw;
        tv{i} = -Rv'*Cv{i};
    end
    varargout{1} = tv;
end
function Xv = to_virtual_camera_coordinate(X,R,C,Rv,Cv)
    if ~isempty(R)
    	Xc = to_camera_coordinate(X,R,C);
    else
        Xc = X;
    end
    Xv = X;
    for i = 1:length(Cv)
        Xv(:,i) = Rv'*Xc(:,i) - Rv'*Cv{i};
    end
end
function p = compute_refractive_point(rs,n,C)
% rs is an array of rays;
% ds is an array of distance (layer);
% n is the common normal
% the following compute the refractive point.
% note: it assumes a common normal for refractive layers
    d = -dot(n,C);
    dots = dot(rs,repmat(n,1,size(rs,2)));% rs is 3xn, n is 3x1
    scales = d ./ dots;% ds is 1xn
    p = C + rs .* scales;
end

