function [X,R,C,timelsq] = refinement_virtual_epipolar(f1s,f2s,n1,n2,d1,d2,ior1,ior2,KK,R0,C0)   
    q1s = imagepoints_to_rays(KK,f1s);
    q2s = imagepoints_to_rays(KK,f2s);
    
    r11 = imagepoints_to_rays(eye(3),double(q1s));
    r12 = compute_refractive_ray(ior1,ior2,r11,n1);
    p1 = compute_refractive_point(r11,d1,n1);
    r21 = imagepoints_to_rays(eye(3),double(q2s));
    r22 = compute_refractive_ray(ior1,ior2,r21,n2);
    p2 = compute_refractive_point(r21,d2,n2);

    % prepare virtual camera parameters
    [Rv1,Cv1] = find_virtual_camera_coordiante(p1,r12,n1);
    [Rv2,Cv2] = find_virtual_camera_coordiante(p2,r22,n2);
            
    pv1 = to_virtual_camera_coordinate(p1,[],[],Rv1,Cv1);
    pv2 = to_virtual_camera_coordinate(p2,[],[],Rv2,Cv2);
    
    % triangulation under refraction
    ray1 = r12; Q1 = p1;
    ray2 = R0 * r22; Q2 = R0 * p2 + C0;
    X = zeros(3,size(f1s,2));
    for ii = 1:size(Q2,2)
        AA = [eye(3) -ray1(1:3,ii) zeros(3,1);eye(3) zeros(3,1) -ray2(1:3,ii)];
        bb = [Q1(1:3,ii);Q2(1:3,ii)];
        sol = AA\bb;
        X(:,ii) = sol(1:3);                
    end   
    
    % because we assume the origin is at the refractive plane.
    R = R0;
    C = C0 + R0*[0;0;d2] - [0;0;d1];    
        
    x0 = [X(:); logSO3(R) ; C];x0 = double(x0);
    
    tic;
    x = solver_gn_own(x0,{Rv1,Rv2},{Cv1,Cv2},{pv1,pv2},{d1,d2},r12,r22,KK);
    timelsq = toc;

    %non-linear optimization
    X = reshape(x(1:end-6),3,[]);    
    R = expSO3(x(end-5:end-3));
    C = x(end-2:end);
end

function varargout = cost_func_reproj(X,R,C,Rv,Cv,pv,fv,omitRt,varargin)
    N = size(X,2);
    M = 3*N + 6;

    tmp1 = R'*X;
    tmp2 = R'*C;
    tmp3 = zeros(3,1);
    if ~isempty(varargin)
        tmp3 = R'*[0;0;-varargin{1}]+[0;0;varargin{2}];%
    end
    
    Xl = tmp1 - tmp2 + tmp3;
    Cvs = horzcat(Cv{:});
    Xv = Rv{1}' * Xl - Rv{1}' * Cvs;
    gv = [fv.*Xv(1,:)./Xv(3,:) - fv.*pv(1,:)./pv(3,:);...
          fv.*Xv(2,:)./Xv(3,:) - fv.*pv(2,:)./pv(3,:)];
    err = gv(:);
    
    fvXv3 = fv./Xv(3,:);
    fvXv1Xv3 = -fv.*Xv(1,:)./(Xv(3,:).^2);
    fvXv2Xv3 = -fv.*Xv(2,:)./(Xv(3,:).^2);
    
    AtA = zeros(M,M);
    Atb = zeros(M,1);
    for i = 1:N
        jac2 = [fvXv3(i) 0 fvXv1Xv3(i);...
                0 fvXv3(i) fvXv2Xv3(i)];
        jac3 = Rv{i}';
        jac4 = zeros(3,M);
        jac4(:,i*3-2:i*3) = R';
        if omitRt == 0
            jac4(:,end-5:end-3) = (skewm(tmp1(:,i)) - skewm(tmp2));%leftJ(-logSO3(R))*
            jac4(:,end-2:end) = -R';
        end
        J = jac2 * jac3 * jac4;
        AtA = AtA + J'*J;
        Atb = Atb + J'*err(i*2-1:i*2);
    end
    varargout{1}=err;
    varargout{2}=AtA;
    varargout{3}=Atb;
end


function varargout = cost_epipolar_own(x,r1,r2,Rv1,Cv1,Rv2,Cv2,d1c,d2c)
    R = expSO3(x(1:3));
    C = x(4:6);
    tt = C-R*[0;0;d2c]+[0;0;d1c];

    N = size(r1,2);
    M = 3*N + 6;    
    
    err = zeros(N,1);

    AtA = zeros(M,M);
    Atb = zeros(M,1);
    J = zeros(1,M);
    for i = 1:size(r1,2)
        rF = R'*skewm((R*Cv2{i}+tt-Cv1{i}));
        err(i) = r2(:,i)'*rF*r1(:,i);
        J(:,end-5:end) = [r2(:,i)'*skewm(rF*r1(:,i))-r2(:,i)'*R'*skewm(r1(:,i))*(-R*skewm(Cv2{i})+R*skewm([0;0;d2c])) -r2(:,i)'*R'*skewm(r1(:,i))];
        AtA = AtA + J'*J;
        Atb = Atb + J'*err(i);
    end
    varargout{1} = err;
    varargout{2}=AtA;
    varargout{3}=Atb;
end

function varargout = cost_func_total_own(x,Rv,Cv,pv,fv,r1,r2,KK)
    X = reshape(x(1:end-6),3,[]);
    R = expSO3(x(end-5:end-3));
    C = x(end-2:end);
    [err1,AtA1,Atb1] = cost_func_reproj(X,eye(3),zeros(3,1),Rv{1},Cv{1},pv{1},KK(1,1),1);
    [err2,AtA2,Atb2] = cost_func_reproj(X,     R,         C,Rv{2},Cv{2},pv{2},KK(1,1),0,fv{1},fv{2});
    [err3,AtA3,Atb3] = cost_epipolar_own(x(end-5:end),r1,r2,Rv{1},Cv{1},Rv{2},Cv{2},fv{1},fv{2});
    AtA = AtA1 + AtA2 + AtA3;
    Atb = Atb1 + Atb2 + Atb3;
    err = [err1;err2;err3];
    varargout{1}=err;
    varargout{2}=AtA;
    varargout{3}=Atb;
end


function x = solver_gn_own(x,Rv,Cv,pv,fv,r1,r2,KK)
    [err,AtA,Atb] = cost_func_total_own(x,Rv,Cv,pv,fv,r1,r2,KK);
   
    errTol = 1e-6;
    found = max(abs(Atb)) < errTol;

    maxIter = 100;
    
    iter = 0;
    epsilon2 = 1e-6;
    tick = 0;
    % routine for gn
%     disp(['iter:',num2str(iter),' cost is:',num2str((err'*err))]);
    while found == false && iter < maxIter
        hlm = (AtA)\(-Atb);
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
            [newerr,newAtA,newAtb] = cost_func_total_own(newx,Rv,Cv,pv,fv,r1,r2,KK);
            err_gain = (err'*err) - (newerr'*newerr);
            
            if abs(err_gain) < epsilon2
                break;
            end
            iter = iter + 1;
%             disp(['iter:',num2str(iter),' cost is:',num2str((newerr'*newerr))]);
            x = newx;
            AtA = newAtA;
            Atb = newAtb;
            found = max(abs(Atb)) < errTol;

            if err_gain<0
                err = newerr;
                tick = tick + 1;
            else
                tick = 0;
            end
            if tick > 3
                break;
            end
        end
    end
end
