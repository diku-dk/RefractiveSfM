function [R,t] = pless_solver_opengv_17pt(imgf1s,imgf2s,KK,n1,n2,d1,d2,ior1,ior2,varargin)
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
    pv1 = to_virtual_camera_coordinate(p1,[],[],Rv1,Cv1);
    
    % prepare virtual camera parameters
    [Rv2,Cv2,~] = find_virtual_camera_coordiante(p2,r22,n2);
    pv2 = to_virtual_camera_coordinate(p2,[],[],Rv2,Cv2);
    
    U1 = pv1(1:2,:)./pv1(3,:);
    U2 = pv2(1:2,:)./pv2(3,:);
    
    K1 = [1 0 0;0 1 0;0 0 1];
    temp = K1 \ [U1;ones(1,size(U1,2))];
    I1_norms = sqrt(sum(temp.*temp));
    I1_normalized = temp ./ repmat(I1_norms,3,1);
    I1_normalized = Rv1{1} * I1_normalized;
    t1 = horzcat(Cv1{:});%zeros(3,size(U,2));    
    
    temp = K1 \ [U2;ones(1,size(U2,2))];
    I2_norms = sqrt(sum(temp.*temp));
    I2_normalized = temp ./ repmat(I2_norms,3,1);
    I2_normalized = Rv2{1} * I2_normalized;
    t2 = horzcat(Cv2{:});%zeros(3,size(U,2));    
    
    [X,~] = opengv('seventeenpt_ransac',[I1_normalized;t1],[I2_normalized;t2]);
%     X = opengv('ge',[I1_normalized;t1],[I2_normalized;t2]);
    
    for i = 1:size(X,3)
        Rs(:,:,i) = X(:,1:3,i);
        ts(:,:,i) = X(:,4,i);
    end
    
    % output the most correct one
    besterr = 1e6;
    bestid = 1;
    for i = 1:size(Rs,3)
        curerr = norm(Rs(:,:,i)'*R1-eye(3),'fro');
        if besterr < curerr
            besterr = curerr;
            bestid = i;
        end
    end
    R = Rs(:,:,bestid);
    t = ts(:,:,bestid);
end
