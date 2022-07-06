function [R,t] = gPnP_Refraction(Q,q,K,n,d,ior1,ior2,varargin)
    r11 = imagepoints_to_rays(K,q);% ideally, should remove distortion, however, rendered images have no distortion, so just do it
    r12 = compute_refractive_ray(ior1,ior2,r11,n);
    p1 = compute_refractive_point(r11,d,n);    

    % prepare virtual camera parameters
    [Rv1,Cv1,tv1] = find_virtual_camera_coordiante(p1,r12,n);
    pv1 = to_virtual_camera_coordinate(p1,[],[],Rv1,Cv1);
    
    U=pv1(1:2,:)./pv1(3,:);
    
    K1 = [1 0 0;0 1 0;0 0 1];
    temp = K1 \ [U;ones(1,size(U,2))];
    I2_norms = sqrt(sum(temp.*temp));
    I2_normalized = temp ./ repmat(I2_norms,3,1);
    t1 = -horzcat(tv1{:});%zeros(3,size(U,2));
%         t1 = horzcat(Cv1{:});%zeros(3,size(U,2));
    X = opengv('gpnp',Q,[I2_normalized;t1]);
    bestid = 0;
    besterr = 1e6;
    for i = 1:size(X,3)
        R2 = X(:,1:3,i);
        t2 = X(:,4,i);
        R1 = R2';
        t1 = -R2'*t2;
        Xreproj = R1 * Q + t1 - t1;
        Xreproj = Xreproj(1:2,:) ./ Xreproj(3,:);
        err = sum(vecnorm(Xreproj - U));
        if err < besterr
            besterr = err;
            bestid = i;
        end
    end
    if bestid == 0
        R = eye(3);
        t = zeros(3,1);
    else
        R2 = X(:,1:3,bestid);
        t2 = X(:,4,bestid);
        R = R2';
        t = -R2'*t2;
    end
    
    Rv = Rv1{1}';
    R  = Rv'*R;
    t = Rv'*t;
end