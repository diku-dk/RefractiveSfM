function [R,t] = gpnp_eccv18_Refraction(Q,q,K,n,d,ior1,ior2,varargin)
%     R1 = varargin{1};
%     tt1 = varargin{2};
    
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

    c1 = horzcat(Cv1{:});
%     I2_normalizedC = Rv1{1}*I2_normalized;

%     Rtruth = Rv1{1}'*R1;
%     ttruth = Rv1{1}'*(-R1*tt1);

    %% self-implementation:
    num_iter = 200;
    logpd = log(1 - 0.9999);
    k = 0;
    max_consensus = 0;
    num_points = 3;
    threshold = 0.5;
    bestR = eye(3);
    bestt = zeros(3,1);
    while k < num_iter
        all_idx = randperm(size(Q,2));
        model_set = all_idx(1:num_points);
        id = model_set;
    
        % set the adapter
        % camera info
        adapter.camera.c2 = t1(:,id(2)) - t1(:,id(1));
        adapter.camera.c3 = t1(:,id(3)) - t1(:,id(1));
        adapter.camera.d1 = I2_normalized(:,id(1));
        adapter.camera.d2 = I2_normalized(:,id(2));
        adapter.camera.d3 = I2_normalized(:,id(3));
        % world info
        adapter.world.p1 = Q(:,id(1));
        adapter.world.p2 = Q(:,id(2));
        adapter.world.p3 = Q(:,id(3));
        [sR,sT,sdepth1,sdepth2,sdepth3] = compute_p3(adapter);
       
        bestid = 0;
        besterr = 1e6;
        N = round(size(sR,2))/3;
        for i = 1:N
            Rsol1 = sR(1:3,3*(i-1)+1:3*i);
            tsol1 = sT(1:3,i);
            if all(imag(Rsol1(:))<1e-6)
                Rsol1 = real(Rsol1);
                tsol1 = real(tsol1);
            else
                continue;
            end
            R1 = Rsol1';
            tsol1 = tsol1 - Rsol1*t1(:,id(1));
            tsol1 = -R1 * tsol1;

            Xreproj = R1 * Q + tsol1 - t1;
            Xreproj = Xreproj(1:2,:) ./ Xreproj(3,:);
            err = sum(vecnorm(Xreproj - U));
            if err < besterr
                besterr = err;
                bestid = i;
            end
        end
        if bestid == 0
            continue;
        else
            Rsol1 = sR(1:3,3*(bestid-1)+1:3*bestid);
            tsol1 = sT(1:3,bestid);
            if all(imag(Rsol1(:))<1e-6)
                Rsol1 = real(Rsol1);
                tsol1 = real(tsol1);
            end
            
            R1 = Rsol1';
            tsol1 = tsol1 - Rsol1*t1(:,id(1));
            tsol1 = -R1 * tsol1;
            Xreproj = R1 * Q + tsol1 - t1;
            Xreproj = Xreproj(1:2,:) ./ Xreproj(3,:);
            model_error = vecnorm(Xreproj - U);
            consensus_idx = find((model_error) < threshold);
            if length(consensus_idx) > max_consensus
                max_consensus = length(consensus_idx);
                % update maxiteration number
                p_inleir = max_consensus / length(model_error);
                new_est_max_iters = round(logpd / (log(1-p_inleir^num_points)+1e-6));
                num_iter = new_est_max_iters;
                
                bestR = R1;
                bestt = tsol1;
            end
        end
        k = k + 1;
    end
        
    Rv = Rv1{1}';
    R  = Rv'*bestR;
    t = Rv'*bestt;
end
