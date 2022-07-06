function [R,t] = refractive_rel_pose_linear_solver(imgf1s,imgf2s,KK,n1,n2,d1,d2,ior1,ior2,varargin)
    q1s = imagepoints_to_rays(KK,imgf1s);
    q2s = imagepoints_to_rays(KK,imgf2s);

    r11 = imagepoints_to_rays(eye(3),double(q1s));
    r12 = compute_refractive_ray(ior1,ior2,r11,n1);
    p1 = compute_refractive_point(r11,d1,n1);
    r21 = imagepoints_to_rays(eye(3),double(q2s));
    r22 = compute_refractive_ray(ior1,ior2,r21,n2);
    p2 = compute_refractive_point(r21,d2,n2);

    if ~isempty(varargin)
        R1 = varargin{1};
    end

    % prepare virtual camera parameters
    [Rv1,Cv1,~] = find_virtual_camera_coordiante(p1,r12,n1);
    
    % prepare virtual camera parameters
    [Rv2,Cv2,~] = find_virtual_camera_coordiante(p2,r22,n2);
    
    t1 = horzcat(Cv1{:});%zeros(3,size(U,2));    
    t2 = horzcat(Cv2{:});%zeros(3,size(U,2));    
    
    A = zeros(size(r12,2),18);
    for i = 1:size(r12,2)
        A(i,:) = [r22(:,i)'*skewm(t2(:,i))*kron(r12(:,i)',eye(3))-r22(:,i)'*kron((skewm(t1(:,i))*r12(:,i))',eye(3)), -r22(:,i)'*kron(r12(:,i)',eye(3))];
    end
    AR = A(:,1:9);
    AE = A(:,10:18);
    ARpinv = pinv(AR);
    N = size(A,1);
    B = (AR*ARpinv - eye(N))*AE;
    [U,S,V] = svd(B);
    sol = V(:,end);
    [mincost, R, t] = solve_kernerl(A,sol);
    
    for i = 1:size(R,3)
         Rt = R(:,1:3,i);
         Rs(:,:,i) = Rt';
%         ts(:,:,i) = X(:,4,i);
    end
%     Rs = R;
%     ts = t;
    
    % output the most correct one
    besterr = 1e6;
    bestid = 1;
    for i = 1:size(Rs,3)
        curerr = norm(Rs(:,:,i)'*R1-eye(3),'fro');
        if curerr < besterr
            besterr = curerr;
            bestid = i;
        end
    end
    
    R = Rs(:,:,bestid);
    
    %% solve t using least square
    A = zeros(size(r12,2),3);
    b = zeros(size(r12,2),1);
    for i = 1:size(r12,2)
        A(i,:) = r22(:,i)'*R'*skewm(r12(:,i));
        b(i,:) = r22(:,i)'*skewm(t2(:,i))*R'*r12(:,i)-r22(:,i)'*R'*skewm(t1(:,i))*r12(:,i);
    end
    t = A\b;
%     t = ts(:,:,bestid);
%     t = -R'*t;
end

function varargout = solve_kernerl(A,sol,varargin)
    sel = 0;
    if ~isempty(varargin)
        sel = 1;
    end
    E = reshape(sol,3,3);
    [Rd,td] = Essential_decomposition(E);
    if sel == 1
        minsumcost = 1e6;
        minid = -1;
        mincost = [];
        for i = 1:size(td,3)
            R0 = Rd(:,:,i);
            t0 = td(:,:,i);
            E0 = skewm(t0)*R0;
            x = [R0(:);E0(:)];
            sumcost = sum(abs(A*x));
            if sumcost < minsumcost
                minsumcost = sumcost;
                minid = i;
                mincost = abs(A*x);
            end
        end
    
        R = Rd(:,:,minid);
        t = td(:,:,minid);
    else
        R = Rd;
        t = td;
        mincost = 0;
    end
    varargout{1} = mincost;
%     varargout{2} = V;
    varargout{2} = R;
    varargout{3} = t;
end