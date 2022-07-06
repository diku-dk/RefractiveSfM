% 8 point estimation algo for planar scene (computes axis, R and vector s)
% find the first two cols of E matrix and s vector using SVD
% completing the last col of E matrix gives 2 E matrices. Both are valid
% thus finally we have 8 rotations and 8 s to check for.
% Find the correction R and vector s using all points



function [bestsol] = EightPointEstimationOneLayerCase1KnownMu(rayDnoise,XYZ,mu,RANSAC_trials,varargin)
    N = size(XYZ,2);

    verbose = 1;
    if ~isempty(varargin)
        verbose = 0;
    end
    
    % get full A matrix
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

        idx = unique(randi(N,8,1));
        while(size(idx,1) < 8)
            idx = unique(randi(N,8,1));
        end

        idx1 = setdiff([1:N]',idx);

        % take 8 points
        K = A(idx,:);

        % find subspace
        [~,~,vv] = svd(K);

        SubSpace = vv(:,9:12);

        % find 4 dimensional subspace for E
        EE = SubSpace(1:9,:);


        % feed it to 5 point algorithm to find coeffs
        coeffs = calibrated_fivepoint_Amit(EE);


        if(isempty(coeffs))
            continue;
        end

        if(isinf(coeffs(1,1)))
            continue;
        end


        % find solutions
        sols = SubSpace*coeffs';

        % find correct scale factor. norm of E is sqrt(2)
        k = size(sols,2);
        for ii = 1:k
            tt = sols(1:9,ii);
            fac = sqrt(2)/norm(tt);
            sols(:,ii) = fac*sols(:,ii);
        end


        % compute error using remaining points
        e = mean(abs(A(idx1,:) * sols));
        e = e(:);
        [~,id] = min(e);
        id = id(1);

        %get correct E value
        sol = sols(:,id);
        clear sols id

        E = reshape(sol(1:9),3,3);
        s = sol(10:12);


        % get normal (axis)
        [uu,~,~] = svd(E);
        ne = uu(:,end);
        if(ne(3) > 0)
            ne = -ne;
        end


        % E will give 4 solutions Compute tz and d for each
        %Compute tz and d

        [Rsol,sSol] = RotationMatricesFromE(E,s,ne);



        %for each of four solutions find tz and d
        d = zeros(4,1);
        t = zeros(3,4);
        for kk = 1:4
            [d(kk,1),t(:,kk)] = EstimateTzOneLayeCase1KnownMuRansac(ne,Rsol(:,:,kk),sSol(:,kk),mu,rayDnoise,XYZ,idx);
        end
        tz = t(3,:)';



        % d should be greater than 0 and tz, translation in z direction should
        % be greater than d
        idx = find(d > 0);
        if(isempty(idx))
            %disp('no solution found for d and tx');
            continue;
        end

        kk = size(idx,1);

        d = d(idx);
        t = t(:,idx);
        Rsol = Rsol(:,:,idx);


        % find the solution with minium 3D error
        C = zeros(kk,1);
        for ii = 1:kk
            C(ii,1) = Compute3DErrorOneLayerCase1(t(:,ii),rayDnoise,XYZ,Rsol(:,:,ii),ne,mu,d(ii,1));
        end
        [val,id] = min(C);



        id = id(1);

        d = d(id);
        t = t(:,id);
        R = Rsol(:,:,id);

        solVec = [ne ; d ; Rot2RPY(R) ; t];

        if(val < BestError)
            BestError = val;
            bestsol = solVec;
            if verbose == 1
                fprintf('3D MSE Error 8 point = %f\n',val);
            end
            if(val < 1e-8)
                break;
            end

        end
    end
end


% estimate tz and depth for one layer given known mu, axis, R and n
% cross t
% case 1 
% do two point ransac

function [d,t] = EstimateTzOneLayeCase1KnownMuRansac(n,R,s,mu,xy,XYZ,idx)
    N = size(XYZ,2);
    % find rotation to make axis along z axis
    w = cross(n,[0;0;-1]);
    w = w/norm(w);
    theta = acos(n'*[0;0;-1]);
    w = theta*w;
    Rc = w2R(w);
    clear w
    clear n

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
        keyboard
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

    idx1 = nchoosek([1:Nsmall],2);

    nn = size(idx1,1);

    bestEr = Inf;
    bestInlier = [];
    for ii = 1:nn
        vecid = idx1(ii,:);
        vecid = vecid(:);
        x = A(vecid,:)\b(vecid);
        d = x(1);
        tz = x(2);
        err = (A*x - b).^2;
        er = median(err);
        if(er < bestEr)
            bestEr = er;
            bestDTz = [d;tz];
            bestInlier = err < 1e-2;
        end
    end
    
    % added, refinement using all inliers
    bestDTz = A(bestInlier,:)\b(bestInlier);

    d = bestDTz(1);
    tz = bestDTz(2);


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















