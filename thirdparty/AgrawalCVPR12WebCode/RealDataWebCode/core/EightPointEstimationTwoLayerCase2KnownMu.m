

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) MERL 2012
% CVPR 2012 Paper Title: A Theory of Multi-Layer Flat Refractive Geometry
% Author: Amit Agrawal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 8 point estimation algorithm along with solving the linear system to find
% calibration parameters
% Assuming known mu and two layer Case 2



function [bestsol] = EightPointEstimationTwoLayerCase2KnownMu(rayDnoise,XYZ,mu,RANSAC_trials)

N = size(XYZ,2);

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


RANSAC_trials

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
        [d(kk,1),t(:,kk)] = EstimateTzTwoLayeCase2KnownMuRansac(ne,Rsol(:,:,kk),sSol(:,kk),mu,rayDnoise,XYZ,idx);
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
    Rsol = Rsol(:,:,idx);
    
    
    % find the solution with minium 3D error
    C = zeros(kk,1);
    for ii = 1:kk
        C(ii,1) = Compute3DErrorTwoLayerCase2(t(:,ii),rayDnoise,XYZ,Rsol(:,:,ii),ne,mu,d(ii,1));
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
        fprintf('3D MSE Error 8 point = %f\n',val);
        
        if(val < 1e-8)
            break;
        end
        
    end
    
    
    
    
    
    
end


