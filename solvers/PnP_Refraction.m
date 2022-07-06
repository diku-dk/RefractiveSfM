function [R,t] = PnP_Refraction(Q,q,K,n,d,ior1,ior2,varargin)
    use_iteropt = 0;
    if ~isempty(varargin)
        use_iteropt = 1;
    end
    r11 = imagepoints_to_rays(K,q);% ideally, should remove distortion, however, rendered images have no distortion, so just do it
    r12 = compute_refractive_ray(ior1,ior2,r11,n);
    p1 = compute_refractive_point(r11,d,n);    

    % prepare virtual camera parameters
    [Rv1,Cv1,tv1] = find_virtual_camera_coordiante(p1,r12,n);
    pv1 = to_virtual_camera_coordinate(p1,[],[],Rv1,Cv1);
    
    U=pv1(1:2,:)./pv1(3,:);
    
    [R,t,err,dims] = REPPnP(Q,U,tv1);

    Rv = Rv1{1}';
    R  = Rv'*R;
    t = Rv'*t;
    
    if use_iteropt == 1
        Rpre = eye(3);
        % refractive points in camera coordinate system, n = [0 0 1]
        norm_ray_u = r11;
        nideal = [0,0,1]';
        if norm(cross(nideal,n))>1e-6
            rotaxis = cross(n,nideal);rotaxis = rotaxis./norm(rotaxis);
            theta = acos(dot(n,nideal));
            Rpre = expSO3(rotaxis.*theta);
            assert(norm(cross(nideal,Rpre*n))<1e-6);
        end
        norm_ray_u_pre = Rpre * norm_ray_u;
        lambda = ior1 / ior2;
        Q1 = [d.*norm_ray_u_pre;norm_ray_u_pre(3,:)];
        R1 = [norm_ray_u_pre(1:2,:).*lambda;sqrt(1-lambda^2+lambda^2.*norm_ray_u_pre(3,:).^2);zeros(1,size(norm_ray_u_pre,2))];
        lparams = zeros(6,size(Q1,2));
        for i = 1:size(Q1,2)
            lparams(:,i) = plucker_line(Q1(:,i),R1(:,i));
        end
        
        %% start nonlinear optimization
        x = [logSO3(R);t];
        P = Q;
        solver = @solver_gn;
        [rvec1,tvec1] = solver(x,P,lparams,Rpre);
        R = expSO3(rvec1);
        t = tvec1;
    end
end

function [R,T,err,varargout] = REPPnP(Pts,impts,tvs)
    sol_iter = 1; %indicates if the initial solution must be optimized
    dims = 4;     %kernel dimensions
    minerror = 0.02; %for the synthetic experiments (f = 800)
    [M, b, Cw, Alph] = PrepareData(Pts,impts,tvs);
    
    Km=robust_kernel_noise(M,b,dims, minerror); %Compute kernel M
    
    dims = 4;
    [R, T, err] = KernelPnP(Cw, Km, dims, sol_iter);
%    T = T - R * mPts;
    varargout{1} = dims;
end

function [R,T,err,varargout] = EPPnP(Pts,impts,tvs)
    sol_iter = 1; %indicates if the initial solution must be optimized
    dims = 4;     %kernel dimensions

    [M, b, Cw, Alph] = PrepareData(Pts,impts,tvs);
    
    Km=kernel_noise(M,b,dims); %Compute kernel M
    
    dims = 4;
    [R, T, err] = KernelPnP(Cw, Km, dims, sol_iter);
%    T = T - R * mPts;
    varargout{1} = dims;
end

function [K, idinliers, i] = robust_kernel_noise(M, b, dimker, minerror)
    m   = size(M,1);
    id  = round(m/8);
    idx = 1:m;

    prev_sv = Inf;
    
    pairs = 1; %each correspondence is a couple of equations

    for i=1:10
        N = M(idx,:);
        bb = b(idx);
        %         [~,~,v] = svd(N'*N);
        if rank(N) == 12
            v = N\bb;
        else
            v = pinv(N)*bb;
        end

        if (pairs)
            error21    = M(1:2:end,:) * v - b(1:2:end);
            error22    = M(2:2:end,:) * v - b(2:2:end);
            error2     = sqrt(error21.^2 + error22.^2);
            
            [sv, tidx] = sort(error2);        

            med = sv(floor(m/8)); 

        else
            error2    = M * v - b;
            [sv, tidx] = sort(error2.^2);
            med = sv(floor(m/2)); 
        end
     
        ninliers = sum(sv<max(med,minerror));

        if (med >= prev_sv)
            break;
        else
            prev_sv = med;
            resv    = N;
            resb    = bb;
            if(pairs)
                residx  = tidx(1:ninliers);
            else
                %always pairs = 1!! :P
               
            end
        end
        
        if(pairs)
            tidx2     = tidx'*2;
            ttidx     = [tidx2-1; tidx2];
            tidx2     = ttidx(:);
            idx       = tidx2(1:2*ninliers);
        else
            idx       = tidx(1:ninliers);
        end
    end
    K=kernel_noise(resv,resb,dimker);
%     K = resv(:,end-dimker+1:end);   
    idinliers = residx;
end

function K=kernel_noise(M,b,dimker)
    K = zeros(size(M,2),dimker);
    [U,S,V] = svd(M);
    K(:,1:dimker-1) = V(:,end-(dimker-2):end);
    if rank(M) < 12
        K(:,end) = pinv(M)*b;
    else
        K(:,end) = pinv(M)*b;
    end
end

function [R, b, mc] = myProcrustes(X,Y)
%X is an structure containing points, points centered in the origin, points
%normalized
%Y are 3D points
    dims = size(Y,2);
    mY = mean(Y,2);
    cY = Y - mY * ones(1,dims);
    ncY = norm(cY(:));
    tcY = cY/ncY;
    
    A = X.nP * tcY';
    [L, D, M] = svd(A);
  
%     R = M * L';
%     
%     if(mY(3)>0 && det(R)<0)
        R = M * diag([1,1,sign(det(M*L'))])* L';
%   end
    
    b = sum(diag(D)) * X.norm/ncY;
    c = X.mP - b*R'*mY;
    mc = c * ones(1,dims);
end

function [R,T, err] = KernelPnP(Cw, Km, dims, sol_iter)

    vK = reshape(Km(:,end),3,dims);
    
    %precomputations
    X.P     = Cw;
    X.mP    = mean(X.P,2);
    X.cP    = X.P - X.mP * ones(1,dims);
    X.norm  = norm(X.cP(:));
    X.nP    = X.cP/X.norm;
    
    %procrustes solution for the first kernel vector
%     if (mean(vK(3,:)<0))
%         vK = -vK;
%     end
    [R,b,mc] = myProcrustes(X,vK);
    
    solV  = b * vK;
    solR  = R;
    solmc = mc;
  
    % procrustes solution using 4 kernel eigenvectors
    if sol_iter
         err = Inf;
         n_iterations=500;
         for iter=1:n_iterations
             % projection of previous solution into the null space
             A = R * (- mc +X.P);
             abcd = Km \ A(:);
             newV = reshape(Km * abcd,3,dims);
             
             %eucliedean error
             newerr = norm(R' * newV + mc - X.P);
%              newerr
             if ((newerr > err) && (iter>2)) || newerr < 1e-6
                 break;
             else
                 %procrustes solution
                 [R,b,mc] = myProcrustes(X,newV);
                 solV = b * newV;
                 
                 solmc = mc;
                 solR = R;
                 err = newerr;
             end
             
         end
    end
       
    R  = solR;
    mV = mean(solV,2);
     
    T = mV - R * X.mP;
end

%% the following are inspired by RePPnP
function [M, b,Cw, Alph] = PrepareData(Pts,impts,tvs,Cw)
    if ~exist('Cw','var')
        Cw=define_control_points()';  
    end
   
    Xw=Pts';
    U=impts;
    
    %compute alphas (linear combination of the control points to represent the 3d points)
    Alph=compute_alphas(Xw,Cw');
   
    %Compute M
    [M,b]=ComputeMb(U(:),Alph,tvs);
end

function Cw=define_control_points()
    % Copyright (C) <2007>  <Francesc Moreno-Noguer, Vincent Lepetit, Pascal Fua>
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
    % Francesc Moreno-Noguer, CVLab-EPFL, September 2007.
    % fmorenoguer@gmail.com, http://cvlab.epfl.ch/~fmoreno/ 
    Cw=[1 0 0;
        0 1 0;
        0 0 1;
        0 0 0];
end

function Alph=compute_alphas(Xw,Cw)
% COMPUTE_ALPHAS Barycentric coordinates computation
%
% Copyright (C) <2007>  <Francesc Moreno-Noguer, Vincent Lepetit, Pascal Fua>
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
% Francesc Moreno-Noguer, CVLab-EPFL, September 2007.
% fmorenoguer@gmail.com, http://cvlab.epfl.ch/~fmoreno/ 
    n=size(Xw,1); %number of 3d points
    %generate auxiliar matrix to compute alphas
    C=[Cw';ones(1,4)];
    X=[Xw';ones(1,n)];
    Alph_=inv(C)*X;
    Alph=Alph_';
end

function [M,b] = ComputeMb(U,Alph,tvs)
    %ATTENTION U must be multiplied by K previously
    M = kron(Alph,[1 0 -1; 0 1 -1]);
    M(:,[[3,6,9,12]]) =  M(:,[3,6,9,12]) .* (U * ones(1,4));
    
    tv = horzcat(tvs{:});
    b1 = [tv(3,:);tv(3,:)];
    b1 = b1(:);
    b1 = b1.*U;
    b2 = [tv(1,:);tv(2,:)];
    b2 = b2(:);
    b = b1 - b2;
end