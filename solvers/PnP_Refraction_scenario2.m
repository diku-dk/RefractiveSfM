function [R,C] = PnP_Refraction_scenario2(Q,q,K,n1w,d1w,ior1,ior2,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code for absolute pose estimation algorithm proposed in:
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
    %% find an idea world coordinate system
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

    Qpre = Rpre * (Q-Cpre);
    n1pre = [0;0;1];
    d1pre = 0;

    % to ideal world coordinate frame.
%     C1pre = Rpre * (C1-Cpre);
%     ref_As_pre = Rpre * (ref_As-Cpre);
%     figure;colormap jet;
%     plot3(C1pre(1),C1pre(2),C1pre(3),'k-d','MarkerSize',15);hold on;
%     plot3(0,0,0,'m-s','MarkerSize',15);hold on;
%     [xx,yy] = meshgrid(min(ref_As_pre(1,:))-1:0.1:max(ref_As_pre(1,:))+1,min(ref_As_pre(2,:))-1:0.1:max(ref_As_pre(2,:))+1);
%     zz = (-d1pre - n1pre(1)*xx - n1pre(2)*yy)/n1pre(3);
%     h3=mesh(xx,yy,zz);hold on
%     h1=plotCameraFrustum(R1*Rpre',-(R1*Rpre')*C1pre,K,0.5,'r',1.5);
%     for ii = 1:num_features
%         h=line([C1pre(1) ; ref_As_pre(1,ii)], [C1pre(2) ; ref_As_pre(2,ii)],[C1pre(3) ; ref_As_pre(3,ii)] ,'LineWidth', 2 , 'Color', 'b');hold on;
%         h.Color(4) = 0.5;
%         h2=plot3(ref_As_pre(1,ii),ref_As_pre(2,ii),ref_As_pre(3,ii),'y-s','MarkerFaceColor','y');hold on;
%         h=line([Qpre(1,ii) ; ref_As_pre(1,ii)], [Qpre(2,ii) ; ref_As_pre(2,ii)],[Qpre(3,ii) ; ref_As_pre(3,ii)] ,'LineWidth', 2, 'Color', 'g');
%         h.Color(4) = 0.5;
%     end
%     h4=plot3(Qpre(1,:),Qpre(2,:),Qpre(3,:),'ko','MarkerFaceColor','k');hold on;
%     xlabel('x: (m)','Interpreter','latex');ylabel('y: (m)','Interpreter','latex');zlabel('z: (m)','Interpreter','latex');
%     title('Refraction in idea world coordinate','Interpreter','latex');
%     legend([h1,h2,h3,h4],{'Camera idea','Refractive Point idea','Refractive Plane idea','World Point idea'},'Interpreter','latex');

    %% initialization with argawal's solutions
    if ~isempty(varargin)
        isplanar = varargin{1};
    else
        isplanar = 0; 
    end
    if isplanar == 1
        [R0,C0] = abs_pose_agw(Q,q,K,ior1,ior2,0,1,0);
    else
        [R0,C0] = abs_pose_agw(Q,q,K,ior1,ior2,0,0,0);
    end
    
    lambda = ior1 / ior2;
    if size(q,1) == 2
        % get homogeneous coordinates
        q = [q;ones(1,size(q,2))];
    end
    % step-1, ray in camera
    ray_u = inv(K)*q;
    norm_ray_u = ray_u ./ sqrt(ray_u(1,:).^2+ray_u(2,:).^2+ray_u(3,:).^2);
    Qpre = Rpre * (Q-Cpre);
    %%%%%%%%%%% above data will not change during optimization %%%%%%%%%%%%
    
    x0 = [logSO3(R0);C0];
    [x,~] = solver_levmar_own(x0,Qpre,norm_ray_u,lambda,Rpre,Cpre);
    R = expSO3(x(1:3));
    C = x(4:6);
    
    
end


function [cost, A, b, J, err] = func(x,Q,norm_ray_u,lambda,Rpre,Cpre)
    Rw = expSO3(x(1:3));% world to camera
    Cw = x(4:6);% camera location in world

    A = zeros(length(x),length(x));
    b = zeros(length(x),1);
    
    % step-4
    ray_w = Rw'*norm_ray_u;% 3xn
    q = Rpre * ray_w;%3xn, ray in the idea world coordinate system
    Ciw = Rpre * (Cw - Cpre);%3x1
    
    if Ciw(3) > 0
        n_iw = [0,0,-1]';
    else
        n_iw = [0,0,1]';
    end
    
    % refractive point in world
    Q1 = [q(3,:).*Ciw + (-Ciw(3))*q;q(3,:)];
    R1 = [q(1:2,:).*lambda;n_iw(3)*sqrt(1-lambda^2+lambda^2.*q(3,:).^2);zeros(1,size(q,2))];
    
    % to plucker
    lparams = zeros(6,size(Q1,2));
    for j = 1:size(Q1,2)
        lparams(:,j) = plucker_line(Q1(:,j),R1(:,j));
    end
    P = Q;
    QcL123 = -cross(P, lparams(1:3,:));
    cost = 0;

    J = zeros(3 * size(Q1,2), 6);
    err = zeros(3 * size(Q1,2), 1);
    for j = 1:size(Q1,2)
        u = QcL123(:,j) + lparams(4:6,j);
        err(j*3-2:j*3) = u;
        cost = cost + u'*u;
        if nargout > 1
            % derivative using chian rule
            dudl123 = -skewm(P(:,j));% done
            AQ = [0 0 0 R1(1,j);0 0 0 R1(2,j);0 0 0 R1(3,j)];% done
            Aq = [-Ciw(3) 0 Ciw(1);0 -Ciw(3) Ciw(2);0 0 0;0 0 1];% done
            
            dqdR = Rpre * skewm(ray_w(:,j));% done
            dQdt = [q(3,j) 0 -q(1,j);0 q(3,j) -q(2,j);0 0 0;0 0 0]*Rpre;% done
            
            AR = [q(3,j) 0 0 0;0 q(3,j) 0 0;0 0 q(3,j) 0];% done
            ARq = [lambda 0 0;0 lambda 0;0 0 n_iw(3)*lambda^2/R1(3,j)*q(3,j);0 0 0];% done
            
            Al1 = dudl123 * ([AQ*Aq*dqdR AQ*dQdt] + ...
                                 [AR*ARq*dqdR zeros(3,3)]);

            
            aa = Q1(1:3,j); bb = R1(1:3,j);
            dl456dQ = [-skewm(bb) zeros(3,1)];
            dl456dR = [skewm(aa) zeros(3,1)];
            
            Al2 = dl456dQ * [Aq*dqdR dQdt] + [dl456dR*ARq*dqdR zeros(3,3)];
            
            % combine
            Atmp = Al1 + Al2;
            A = A + Atmp'*Atmp;
            b = b + Atmp'*u;
            J(j*3-2:j*3,:) = Atmp;
        end
    end        
end

function [x,finalcost] = solver_levmar_own(x,Q,norm_ray_u,lambda,Rpre,Cpre)
    k = 0;
    v = 2;
    
    J = [];
    err = [];
    
    [~,~,~,J,err] = func(x,Q,norm_ray_u,lambda,Rpre,Cpre);
    A = (J'*J); g = (J'*err);
    errTol = 1e-10;
    found = max(abs(g)) < errTol;
    tau = 1e-3;
    miu = tau * max(diag(A));
    maxIter = 1000;
    iter = 0;
    epsilon2 = 1e-10;
     % routine for lm
     newerr = err;
    disp(['init cost is:',num2str((newerr'*newerr))]);
%     return;
    while found == false && iter < maxIter
        iter = iter + 1;
        hlm = (A+miu*eye(size(A,1),size(A,2)))\(-g);
        hnorm = sqrt(hlm'*hlm);
        if hnorm <= epsilon2*(sqrt(x'*x)+epsilon2)
            found = true;
        else
            R = expSO3(x(end-5:end-3));
            C = x(end-2:end);
            newR = expSO3(logSO3(R))*expSO3(hlm(end-5:end-3));
            newC = C + hlm(end-2:end);
            newx = [logSO3(newR) ; newC];
            [~,~,~,newJ,newerr] = func(newx,Q,norm_ray_u,lambda,Rpre,Cpre);
            err_gain = (err'*err) - (newerr'*newerr);
            gain = err_gain / (0.5*(hlm'*(miu*hlm-g)));
            if gain > 0
                disp(['iter:',num2str(iter),' cost is:',num2str((newerr'*newerr))]);
                x = newx;
                err = newerr;
                A = (newJ'*newJ);
                g = (newJ'*newerr);
%                 A = newA;
%                 g = newg;
                found = max(abs(err)) < errTol;
                miu = miu * max(1/3,1-(2*gain-1)^3);v = 2;
                oldcost = sqrt(newerr'*newerr);
            else
                miu = miu * v;v = v*2;
            end
        end
    end
    disp(['final cost is:',num2str((newerr'*newerr))]);
    if nargout > 1
        finalcost=sqrt(newerr'*newerr);
    end
end

function [R,C] = abs_pose_agw(Q,q,K,ior1,ior2,varargin)
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






