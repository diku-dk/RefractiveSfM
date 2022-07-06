clear all;close all;clc;

addpath(genpath('../'));
mu = 1.5 ; % refractive index of medium
COP = [0;0;0];  % Camera coordinate system is origin

noise_level = [0 0.1 0.5 1 2 3 4 5];
sim_num = 100;
KK = [4800 0 960;0 4800 540;0 0 1];
feature_num = 100;

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);


if ~exist('./data/res.mat','file')
    fbar = waitbar(0,'Please wait...');
    warning off;
    for i = 1:length(noise_level)
        stdnoise = noise_level(i);
        for j = 1:sim_num
            % generate a flat refractive layer at distance d0
            d0 = 10;
            % generate random normal of the flat refractive layer
            n = [randn(1)*0.3 ; randn(1)*0.3 ; -1 + randn(1)*0.3];
            n = n/sqrt(n'*n);
            if(n(3)>=0)
                n = -n;
            end
            
            waitbar(((i-1)*sim_num+j)/(length(noise_level)*sim_num),fbar,'Processing...');

            % Now we generate image points and find the final refracted ray working
            % from the camera.
            deltad = 1+rand(1)*10;

            for trial = 1:feature_num
                % generate random image feature point in normalized coordinate system
                x2D = rand(2,1) - 0.5;
                rayD = [x2D ; 1];
                %find intersection of this ray with the flat refracting plane
                lamda = -d0/(n'*rayD);
                % lamda should be > 0 so that 3D point is in front
                while(lamda < 0)
                    x2D = rand(2,1) - 0.5;
                    rayD = [x2D ; 1];
                    %find intersection of this ray with the flat plane
                    % lamda*vi should lie on the plane
                    lamda = -d0/(n'*rayD);
                end
                q = lamda*rayD; % q lies on flat refracting plane
                % find refracted ray
                [v2,~,~] = RefractedRay(q,n,1,mu);
                %find a 3D point on this ray at some random distance
                fac = abs((d0+deltad-(n'*q))/(n'*v2)*norm(v2));

                p = q + fac*v2/norm(v2);

                x2I = KK * rayD;

                xy(:,trial) = x2I(1:2) + randn(2,1)*stdnoise;
                XYZ(:,trial) = p;
            end

            Rgt = eye(3);
            Tgt = zeros(3,1);
            rotaxis = cross(n,[0,0,1]');rotaxis = rotaxis./norm(rotaxis);
            rotangle = acos(dot(n,[0,0,1]'));
            if abs(rotangle) > 1e-6
                Rgt = expSO3(rotangle.*rotaxis)';
                XYZnew = Rgt' * XYZ;
                Tgt = -Rgt*[0;0;-mean(XYZnew(3,:))];
                XYZ = XYZnew+[0;0;-mean(XYZnew(3,:))];
            end
            w = Rot2RPY(Rgt);

            disp('============1===========')
            CameraRays = [xy;ones(1,size(xy,2))];
            Pall = XYZ;
            N = size(Pall,2);
            ray2D = inv(KK)*CameraRays;
            try
            tic;
            solVec = AxisEstimationPlanarSceneOneLayerCase1(ray2D,Pall,mu);
            time8pt = toc;
            catch
                j = j - 1;
                continue;
            end
            ne8pt = solVec(1:3);
            tau8pt = solVec(4);
            R8pt = RPY2Rot(solVec(5:7));
            t8pt = solVec(8:10);
            w8pt = solVec(5:7);

            errne8pt(i,j) = acos(dot(ne8pt,n)./norm(ne8pt)./norm(n));
            errtau8pt(i,j) = abs(d0 - tau8pt);
            errR8pt(i,j) = norm(Rgt'*R8pt-eye(3),'fro');
            errt8pt(i,j) = norm(Tgt - t8pt);
            runtime8pt(i,j) = time8pt;

            disp('============2===========')
            %non-linear optimization
            x0 = [ne8pt(1:2)/abs(ne8pt(3)) ; w8pt ; t8pt ; tau8pt ; ];
            options = optimset('Display','off','TolFun',1e-6,'TolX',1e-6,'MaxIter',5000,'MaxFunEvals',10000);
            tic;
            x = lsqnonlin(@CostFuncOneLayersCase1,x0,[],[],options,CameraRays,Pall,KK,mu);
            timelsq = toc;
            ne8pt_1 = [x(1:2) ; -1];
            ne8pt_1 = ne8pt_1/norm(ne8pt_1);
            w8pt_1 = x(3:5);
            t8pt_1 = x(6:8);
            tau8pt_1 = x(9);
            R8pt_1 = RPY2Rot(w8pt_1);
            errnelsq(i,j) = acos(dot(ne8pt_1,n)./norm(ne8pt_1)./norm(n));
            errtaulsq(i,j) = abs(d0 - tau8pt_1);
            errRlsq(i,j) = norm(Rgt'*R8pt_1-eye(3),'fro');
            errtlsq(i,j) = norm(Tgt - t8pt_1);
            runtimelsq(i,j) = timelsq;

            disp('==============3=================')
            refra_n1 = 1;
            refra_n2 = mu;
            R1 = RPY2Rot(solVec(5:7));
            t1 = -R1'*solVec(8:10);
            views = 1;
            % covert normal to Rpre
            n1t = R1'*ne8pt;n1t = n1t./norm(n1t);
            d = abs(dot(n1t,t1)) - tau8pt;
            % covert normal to Rpre
            % n1t = R1'*ne8pt;n1t = n1t./norm(n1t);
            if n1t(3)<0, n1t = n1t*-1;end
            rotaxis = cross(n1t,[0;0;1]);
            rotaxis = rotaxis ./ norm(rotaxis);
            theta = acos(dot(n1t,[0;0;1]));
            if abs(theta) < 1e-6
                Rpre = eye(3);
            else
                Rpre = expSO3(theta*rotaxis);
            end
            assert(abs(dot(Rpre*n1t,[0;0;1])-1)<1e-6);
            tic;
            [nopt,dopt,R,t] = solveJointOptimization({XYZ},{xy}, KK, refra_n1, refra_n2, Rpre, d, {R1}, {t1},views);
            timeopt = toc;
            x(3:5) = Rot2RPY(R);
            x(6:8) = t;
            tmp = -R'*t;
            % x(9) = tmp(3) - dopt;
            x(9) = abs(dot(nopt,tmp)) - dopt;
            ne8pt_2 = R*nopt;
            w8pt_2 = x(3:5);
            t8pt_2 = x(6:8);
            tau8pt_2 = x(9);

            errneopt(i,j) = acos(dot(ne8pt_2,n)./norm(ne8pt_2)./norm(n));
            errtauopt(i,j) = abs(d0 - tau8pt_2);
            errRopt(i,j) = norm(Rgt'*R-eye(3),'fro');
            errtopt(i,j) = norm(Tgt - t8pt_2);
            runtimeopt(i,j) = timeopt;
        end
    end
    waitbar(1,fbar,'End...');
    close(fbar);
    save('./data/res.mat','errne8pt','errtau8pt','errR8pt','errt8pt','runtime8pt',...
        'errnelsq','errtaulsq','errRlsq','errtlsq','runtimelsq',...
        'errneopt','errtauopt','errRopt','errtopt','runtimeopt');
else
    load('./data/res.mat');
    
    res_error_n = cell(length(noise_level),3);
    res_error_d = cell(length(noise_level),3);
    res_error_r = cell(length(noise_level),3);
    res_error_t = cell(length(noise_level),3);
    
    for j = 1:length(noise_level)
        res_error_r{j,  1} = errR8pt(j,:);
        res_error_t{j,  1} = errt8pt(j,:);
        res_error_n{j,  1} = errne8pt(j,:);
        res_error_d{j,  1} = errtau8pt(j,:);
        
        res_error_r{j,  2} = errRlsq(j,:);
        res_error_t{j,  2} = errtlsq(j,:);
        res_error_n{j,  2} = errnelsq(j,:);
        res_error_d{j,  2} = errtaulsq(j,:);
        
        res_error_r{j,  3} = errRopt(j,:);
        res_error_t{j,  3} = errtopt(j,:);
        res_error_n{j,  3} = errneopt(j,:);
        res_error_d{j,  3} = errneopt(j,:);
    end
    xlabels = categorical(convertStringsToChars(string(noise_level)));
    box_labels = convertStringsToChars(string(noise_level));
    plot_case = {'8pt','forward','backward'};
    cmap = lines(3);
    cmapalpha = [cmap 0.3*ones(size(cmap,1),1)];
    figure
    multiple_boxplot(res_error_r, box_labels, plot_case, cmapalpha');
    title('Error of R: $\|R_{th}^T\hat{R}\|_{fro}$','Interpreter','latex');
    xlabel('Noise Level');
    ylabel('Error');
    figure
    multiple_boxplot(res_error_t, box_labels, plot_case, cmapalpha');
    title('Error of t: $\|t_{th}-\hat{t}\|$','Interpreter','latex');
    xlabel('Noise Level');
    ylabel('Error');
    figure
    multiple_boxplot(res_error_n, box_labels, plot_case, cmapalpha');
    title('Error of normal: $\arccos(\frac{\|n_{th}\cdot \hat{n}\|}{\|n_{th}\| \|\hat{n}\|})$','Interpreter','latex');
    xlabel('Noise Level');
    ylabel('Error');
    figure
    multiple_boxplot(res_error_d, box_labels, plot_case, cmapalpha');
    title('Error of distance: $\|d_{th}-\hat{d}\|$','Interpreter','latex');
    xlabel('Noise Level');
    ylabel('Error');
    
    figure
    plot(noise_level,mean(runtime8pt,2),'-o','Color',cmapalpha(1,:));
    plot(noise_level,mean(runtimelsq,2),'-s','Color',cmapalpha(2,:));
    plot(noise_level,mean(runtimeopt,2),'-d','Color',cmapalpha(3,:));
    title('Runtime');
end


function [n,d,R,t] = solveJointOptimization(Q,q,K,ior1,ior2,Rpre,d,R0,t0,id)
% this n is norm of the refractive plane in camera coordinate system
% d is the perpendicular distance between the refractive plane to the
% camera center.

    %% here starts joint optimization
    views = length(Q);
    views = min(views,id);
    rvec=cell(views,1);
    tvec=cell(views,1);
    norm_ray_u = cell(views,1);
    R = cell(views,1);
    t = cell(views,1);
    lambda = ior1 / ior2;

    for i = 1:views
        % get initial pose
%         [rvec{i},tvec{i}] = cv.solvePnP(Q{i}',q{i}',K,'Method','Iterative');
        R{i} = R0{i};
        t{i} = t0{i};% t in world coordinate
        if size(q{i},1) == 2
            % get homogeneous coordinates
            q{i} = [q{i};ones(1,size(q{i},2))];
        end
        % step-1, ray in camera
        ray_u = inv(K)*q{i};
        norm_ray_u{i} = ray_u ./ sqrt(ray_u(1,:).^2+ray_u(2,:).^2+ray_u(3,:).^2);
    end
    %%%%%%%%%%% above data will not change during optimization %%%%%%%%%%%%

    len = views*6+3+1;
    x0 = zeros(len,1);
    for i = 1:views
        x0(i*6-5:i*6) = [logSO3(R{i});t{i}];
    end
    x0(len-3:end) = [logSO3(Rpre);d];

    [x,~] = solver_gn(x0,Q,norm_ray_u,lambda,views);

    Ropt = cell(views,1);
    topt = cell(views,1);
    for i = 1:views
        Ropt{i} = expSO3(x(i*6-5:i*6-3));
        topt{i} = x(i*6-2:i*6);
    end

    R = Ropt{1};
    t = -R*topt{i};

    Rpre = expSO3(x(end-3:end-1));
    d = x(end);
    % refractive points in camera coordinate system, n = [0 0 1]
    nideal = [0,0,1]';
    n = Rpre'*nideal;
end

function [costs, A, b] = func(x,Q,norm_ray_u,lambda,views)
    Rw = cell(views,1);
    tw = cell(views,1);
    for i = 1:views
        Rw{i} = expSO3(x(i*6-5:i*6-3));
        tw{i} = x(i*6-2:i*6);
    end
    Rpre = expSO3(x(end-3:end-1));
    d = x(end);

    costs = 0;
    A = zeros(length(x),length(x));
    b = zeros(length(x),1);
    for i = 1:views
        ray_w = Rw{i}'*norm_ray_u{i};% 3xn
        q = Rpre * ray_w;%3xn
        t = Rpre * tw{i};%3x1
        % refractive point in world
        Q1 = [q(3,:).*t + (-t(3)+d)*q;q(3,:)];
        R1 = [q(1:2,:).*lambda;-sqrt(1-lambda^2+lambda^2.*q(3,:).^2);zeros(1,size(q,2))];
        % to plucker
        lparams = zeros(6,size(Q1,2));
        for j = 1:size(Q1,2)
            lparams(:,j) = plucker_line(Q1(:,j),R1(:,j));
        end

        P0 = Q{i};
        P = Rpre * P0;
        QcL123 = -cross(P, lparams(1:3,:));

        % single view cost and grad
        cost = 0;
        A0 = zeros(10,10);
        b0 = zeros(10,1);
        for j = 1:size(Q1,2)
            u = QcL123(:,j) + lparams(4:6,j);
            cost = cost + u'*u;

            if nargout > 1
                % derivative using chian rule
                dudp = skewm(lparams(1:3,j));
                dpdRpre = -Rpre*skewm(P0(:,j));
                Ap = dudp * [zeros(3,3) zeros(3,3) dpdRpre zeros(3,1)];


                dudl123 = -skewm(P(:,j));
                AQ = [0 0 0 R1(1,j);0 0 0 R1(2,j);0 0 0 R1(3,j)];
                Aq = [-t(3)+d 0 t(1);0 -t(3)+d t(2);0 0 d;0 0 1];
                dqdR = Rpre * skewm(ray_w(:,j));
                dqdRpre = -Rpre * skewm(ray_w(:,j));
                dQdt = [q(3,j) 0 -q(1,j);0 q(3,j) -q(2,j);0 0 0;0 0 0];
                dQdd = [q(1,j);q(2,j);q(3,j);0];

                AR = [q(3,j) 0 0 0;0 q(3,j) 0 0;0 0 q(3,j) 0];
                ARq = [lambda 0 0;0 lambda 0;0 0 -lambda^2/R1(3,j)*q(3,j);0 0 0];
                Al1 = dudl123 * ([AQ*Aq*dqdR AQ*dQdt AQ*Aq*dqdRpre AQ*dQdd] + ...
                                 [AR*ARq*dqdR zeros(3,3) AR*ARq*dqdRpre zeros(3,1)]);

                aa = Q1(1:3,j); bb = R1(1:3,j);
                dl456dQ = [-skewm(bb) zeros(3,1)];
                dl456dR = [skewm(aa) zeros(3,1)];


                Al2 = dl456dQ * [Aq*dqdR dQdt Aq*dqdRpre dQdd] + [dl456dR*ARq*dqdR zeros(3,3) dl456dR*ARq*dqdRpre zeros(3,1)];

                % combine
                Atmp = Ap + Al1 + Al2;
                A0 = A0 + Atmp'*Atmp;
                b0 = b0 + Atmp'*u;
            end
        end
        costs = costs + cost;
        A(i*6-5:i*6,i*6-5:i*6) = A(i*6-5:i*6,i*6-5:i*6) + A0(1:6,1:6);
        A(i*6-5:i*6,end-3:end) = A(i*6-5:i*6,end-3:end) + A0(1:6,7:10);
        A(end-3:end,i*6-5:i*6) = A(end-3:end,i*6-5:i*6) + A0(7:10,1:6);
        A(end-3:end,end-3:end) = A(end-3:end,end-3:end) + A0(7:10,7:10);
        b(i*6-5:i*6) = b(i*6-5:i*6) + b0(1:6);
        b(end-3:end) = b(end-3:end) + b0(7:10);
    end
end

function [x,finalcost] = solver_gn(x,Q,norm_ray_u,lambda,views)
    [oldcost, A, b] = func(x,Q,norm_ray_u,lambda,views);

    maxiter = 1000;
    exitcnt = 0;
    for iter = 1:maxiter
        % call cost
        % compute gradient
%         disp(['iter', num2str(iter), ', cost: ', num2str(oldcost)]);
        xi = -(A)\(b);

        newx = x;
        for i = 1:views
            newx(i*6-5:i*6-3) = logSO3(expSO3(x(i*6-5:i*6-3))*expSO3(xi(i*6-5:i*6-3)));
            newx(i*6-2:i*6) = x(i*6-2:i*6) + xi(i*6-2:i*6);
        end
        newx(end-3:end-1) = logSO3(expSO3(x(end-3:end-1))*expSO3(xi(end-3:end-1)));
        newx(end) = newx(end) + xi(end);

        [newcost, newA, newb] = func(newx,Q,norm_ray_u,lambda,views);
        if abs(newcost - oldcost) < 1e-10
            break;
        end
        if newcost < oldcost
            oldcost = newcost;
            x = newx;
            A = newA;
            b = newb;
            exitcnt = 0;
        else
%             oldcost = newcost;
            x = newx;
            A = newA;
            b = newb;
            exitcnt = exitcnt + 1;
            if exitcnt > 10
                break;
            end
        end
    end

    if nargout > 1
        finalcost=oldcost;
    end
end





