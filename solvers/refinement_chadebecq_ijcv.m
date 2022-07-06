%%%%%%%%%%%%%%%%%%%%%% end check optimization performance %%%%%%%%%%%%%%%
function [X,R,C,timelsq] = refinement_chadebecq_ijcv(f1s,f2s,n1,n2,d1,d2,ior1,ior2,miu,KK,R0,C0)        
    q1s = imagepoints_to_rays(KK,f1s);
    q2s = imagepoints_to_rays(KK,f2s);
    
    r11 = imagepoints_to_rays(eye(3),double(q1s));
    r12 = compute_refractive_ray(ior1,ior2,r11,n1);
    p1 = compute_refractive_point(r11,d1,n1);
    r21 = imagepoints_to_rays(eye(3),double(q2s));
    r22 = compute_refractive_ray(ior1,ior2,r21,n2);
    p2 = compute_refractive_point(r21,d2,n2);
    
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

    % data has been prepared for chadebecq ijcv soluion.
    [NrP, Chats] = prepare_data_for_chadebecq(1/miu,d1,d2,size(X,2),q2s);

    % because we assume the origin is at the refractive plane.
    R = R0;
    C = C0 + R0*[0;0;d2] - [0;0;d1];    
    
    %non-linear optimization
    x0 = [X(:); logSO3(R) ; C];x0 = double(x0);
    options = optimset('Display','off','TolFun',1e-6,'TolX',1e-6,'MaxIter',1000,'MaxFunEvals',10000,'Algorithm','levenberg-marquardt');
    tic;
    x = lsqnonlin(@cost_func_total,x0,[],[],options,f1s,f2s,n1,n2,d1,d2,q1s,q2s,NrP,Chats,miu,KK);
    timelsq = toc;
    X = reshape(x(1:end-6),3,[]);    
    R = expSO3(x(end-5:end-3));
    C = x(end-2:end);
end


function [err] = cost_epipolar(x, q1,q2,NrP,Chats)
    R = expSO3(x(1:3));
    C = x(4:6);

    Ss = zeros(21,36);
    ind = [1,1;2,2;2,7;3,8;4 3;4 13;5 9;5 14;6 15;7 4;7 19;8 10;8 20;9 16;9 21;10 22;11 5;11 25;12 11;12 26;13 17;13 27;14 23;14 28;15 29;16 6;16 31;17 12;17 32;18 18;18 33;19 24;19 34;20 30;20 35;21 36];
    Ss(ind(:,1)+(ind(:,2)-1)*21) = 1;
    Ds=(diag(sum(Ss,2)));

    T = [R zeros(3,3);skewm(C)*R R];
    T2 = [T(6,[6,1,2,4,5,3]);...
          T(1,[6,1,2,4,5,3]);...
          T(2,[6,1,2,4,5,3]);...
          T(4,[6,1,2,4,5,3]);...
          T(5,[6,1,2,4,5,3]);...
          T(3,[6,1,2,4,5,3])];
    T2hat = inv(Ds)*Ss*kron(T2,T2)*Ss';
    
    T2NrP = T2hat'*NrP;
    
    fhat6 = @(q) ([q(1)^2 q(1)*q(2) q(2)^2 q(1)*q(3) q(2)*q(3) q(3)^2 q(1)*q(4) q(2)*q(4) q(3)*q(4) q(4)^2 q(1)*q(5) q(2)*q(5) q(3)*q(5) q(4)*q(5) q(5)^2 q(1)*q(6) q(2)*q(6) q(3)*q(6) q(4)*q(6) q(5)*q(6) q(6)^2]');
    fhat = @(q) ([q(1)*q(1);q(1)*q(2);q(2)*q(2);q(1)*q(3);q(2)*q(3);q(3)*q(3)]);
    
    N = size(q1,2);
    err = zeros(N,1);
    for i = 1:N
        Chat = Chats{i};
        rF = Chat'*T2NrP;
        err(i) = fhat6([q2(:,i);q2(:,i)/q2(3,i)])' * rF * [fhat(q1(1:3,i))/(q1(3,i)^2);fhat(q1(1:3,i))];
    end
end

function cost = cost_func_reproj(R,C,n,d,features,X,KK,miu,d1,d2)
    Cnew = C;
    if ~isempty(d1)
        Cnew = C + [0;0;d1] - R*[0;0;d2];
    end
    cost = CostFuncOneLayersCase1(R,Cnew,n,d,features,X,KK,miu);
end

function [err] = cost_func_total(x,f1s,f2s,n1,n2,d1,d2,q1,q2,NrP,Chats,miu,KK)
    X = reshape(x(1:end-6),3,[]);
    R = expSO3(x(end-5:end-3));
    C = x(end-2:end);
   
    err1 = cost_func_reproj(eye(3),zeros(3,1),-n1,d1,f1s,X,KK,miu,[],[]);
    err2 = cost_func_reproj(R,C,-n2,d2,f2s,X,KK,miu,d1,d2);
    err3 = cost_epipolar(x(end-5:end),q1,q2,NrP,Chats);
    
    err = [err1;err2;err3];
end