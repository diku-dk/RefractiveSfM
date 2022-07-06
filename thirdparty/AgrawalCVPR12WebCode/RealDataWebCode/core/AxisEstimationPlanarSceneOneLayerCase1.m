% 8 point estimation algo for planar scene (computes axis, R and vector s)
% find the first two cols of E matrix and s vector using SVD
% completing the last col of E matrix gives 2 E matrices. Both are valid
% thus finally we have 8 rotations and 8 s to check for.
% Find the correction R and vector s using all points


function [bestsol] = AxisEstimationPlanarSceneOneLayerCase1(rayDnoise,XYZ,mu)

    N = size(XYZ,2);

    A = zeros(N,12);
    for ii = 1:N
        P = XYZ(:,ii);
        v = rayDnoise(:,ii);
        A(ii,:) = [ kron(P',v') v'];
    end

    BestError = Inf;

    if(~exist('RANSAC_trials','var'))
        RANSAC_trials = 300;
    end


    RANSAC_trials;

    for trial = 1:RANSAC_trials

        % take 8 points
        idx = unique(randi(N,8,1));
        while(size(idx,1) < 8)
            idx = unique(randi(N,8,1));
        end

        Asmall = A(idx,[1:6 10:12]);

        % find svd
        [~,~,vv] = svd(Asmall,0);
        sol = vv(:,end);

        % get 2 cols of E matrix
        s = sol(7:9);
        x = reshape(sol(1:6),3,2);
        clear sol

        % find E matrices. get two real solutions differing in sign
        Ef = SolveLastColOfEmatrix(x);

        if(isempty(Ef))
            continue;
        end

        kk = size(Ef,3);

        for ii = 1:kk
            tt = Ef(:,:,ii);
            tt = tt(:);
            fac = sqrt(2)/norm(tt);
            Ef(:,:,ii) = fac*Ef(:,:,ii);
        end
        s = fac*s;

        % each E matrix will give 4 rotations
        AllR = zeros(3,3,4*kk);
        Alls = zeros(3,4*kk);

        for ii = 1:kk
            [uu,~,~] = svd(Ef(:,:,ii));
            ne = uu(:,end);
            if(ne(3) > 0)
                ne = -ne;
            end
            [RR,ss] = RotationMatricesFromE(Ef(:,:,ii),s,ne);
            AllR(:,:,(ii-1)*4+1:ii*4) = RR;
            Alls(:,(ii-1)*4+1:ii*4) = ss;
        end

        %for each solution find tz and d
        nn = size(Alls,2);
        d = zeros(nn,1);
        t = zeros(3,nn);
        for kk = 1:nn
            [d(kk,1),t(:,kk)] = EstimateTzOneLayeCase1KnownMu(ne,AllR(:,:,kk),Alls(:,kk),mu,rayDnoise,XYZ,idx);
        end

        tz = t(3,:)';

        % d should be greater than 0 and tz, translation in z direction should
        % be greater than d
        idx = find(d > 0 & tz > d);
        if(isempty(idx))
            %disp('no solution found for d and tx');
            continue;
        end

        kk = size(idx,1);



        d = d(idx);
        t = t(:,idx);
        Rsol = AllR(:,:,idx);

        % find solution with minimum error using all points
        C = zeros(kk,1);
        for ii = 1:kk
            C(ii,1) = Compute3DErrorOneLayerCase1(t(:,ii),rayDnoise,XYZ,Rsol(:,:,ii),ne,mu,d(ii,1));
        end
        [val,id] = sort(C);


        % check smallest two. We will always get two solutions for R which has
        % same error, since last column of R does not matter in R*P + t
        % then check which R matrix is valid rotation matrix (other will be a
        % reflection)
        R1 = Rsol(:,:,id(1));
        R2 = Rsol(:,:,id(2));
        e1 = sum(abs(cross(R1(:,1),R1(:,2)) - R1(:,3)));
        e2 = sum(abs(cross(R2(:,1),R2(:,2)) - R2(:,3)));


        if(e1 < e2)
            id = id(1);
            val = val(1);
        else
            id = id(2);
            val = val(2);
        end

        d = d(id);
        t = t(:,id);
        R = Rsol(:,:,id);

        solVec = [ne ; d ; Rot2RPY(R) ; t];

        if(val < BestError)
            BestError = val;
            bestsol = solVec;
%             fprintf('3D MSE Error 8 point planar scene  = %f\n',val);

            if(val < 1e-6)
                break;
            end

        end

    end
end




% estimate tz and depth for a one layer given known mu, axis, R and vector s
% case 1
function [d,t] = EstimateTzOneLayeCase1KnownMu(n,R,s,mu,xy,XYZ,idx)

    N = size(XYZ,2);

    % find rotation to make axis along -z axis
    w = cross(n,[0;0;-1]);
    w = w/norm(w);
    theta = acos(n'*[0;0;-1]);
    w = theta*w;
    Rc = w2R(w);
    clear w

    clear n

    %check Rc*n should be [0;0;-1]


    % rotate s = cross(n,t)
    ss = Rc * s;
    % now tx and tz can be obtained. We need to find tz
    t = [-ss(2) ; ss(1) ; 0];
    clear ss


    % rotate 3D points
    XYZ = Rc*R*XYZ + t*ones(1,N);

    n2D = [0;-1];

    % find z1, z2
    z1 = [0;0;1];

    %
    zp = [0;1];

    if(isempty(idx))
        idx = [1:N]';
    end


    Nsmall = size(idx,1);

    A = zeros(Nsmall,2);
    b = zeros(Nsmall,1);

    for jj = 1:Nsmall

        ii = idx(jj,1);

        % camera ray
        v = Rc*xy(:,ii);

        P = XYZ(:,ii);

        %find z2
        z2 = cross(z1,cross(z1,v));
        z2 = z2/norm(z2);
        % make along x axis
        if(z2(1)<0)
            z2 = -z2;
        end

        %find projection
        vp0 = [z2'*v ; z1' * v];
        vp0 = vp0/norm(vp0);

        %find refracted ray
        [vp1,~,~] = RefractedRay(vp0,n2D,1,mu);

        u = [z2'*P ; z1'*P];

        %make linear system
        b(jj,1) = -Cross2(vp1,u);

        A(jj,:) = [Cross2(vp1,vp0)/(vp0'*n2D) Cross2(vp1,zp)];

    end

    x = A\b;

    d = x(1);
    tz = x(2);

    % put tz in t
    t(3) = tz;

    % rotate back into original coordinate system
    t = Rc'*t;
end



% Compute 3D error between the back-projected ray and the given 3D point
% for one layers Case 1
function [C] = Compute3DErrorOneLayerCase1(t,xy,XYZ,R,n,mu,d1)

    N = size(XYZ,2);

    % transform points
    P = R*XYZ + t*ones(1,N);

    e = zeros(N,1);

    for ii = 1:N
        pt = P(:,ii);
        v = xy(:,ii);

        q1 = -d1*v/(v'*n);


        [v1,~,~] = RefractedRay(v,n,1,mu);


        %for case 1, refracted ray is v1, q1-pt should have the same direction as v1

        e(ii,1) = norm(cross(v1,q1 - pt))/norm(v);

    end

    C = median(e.^2);
end















