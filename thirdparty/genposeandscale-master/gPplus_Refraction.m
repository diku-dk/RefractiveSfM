function [R,t] = gPplus_Refraction(Q,q,K,n,d,ior1,ior2,varargin)
%     R1 = varargin{1};
%     t1 = varargin{2};
    
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
    t1 = -horzcat(tv1{:});

%     Rtruth = Rv1{1}'*R1;
%     ttruth = Rv1{1}'*(-R1*t1);
    
    [Psolns,ssolns] = genposeandscale(t1,I2_normalized,Q);
    
    bestid = 0;
    besterr = 1e6;
    for i = 1:length(Psolns)
        Rsol1 = Psolns{i}(1:3,1:3);
        Rsol1 = Rsol1';
        tsol1 = Psolns{i}(1:3,4);
        tsol1 = -Rsol1 * tsol1;
        ssol1 = ssolns(i);
        
        Xreproj = ssol1*(Rsol1 * Q + tsol1) - t1;
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
        R2 = Psolns{bestid}(1:3,1:3);
        t2 = Psolns{bestid}(1:3,4);
        ssol1 = ssolns(bestid);
        R = R2';
        t = -R * t2;
        t = t./ssol1;
    end
    
    Rv = Rv1{1}';
    R  = Rv'*R;
    t = Rv'*t;
end
