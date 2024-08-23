function [R,t] = Essential_decomposition(varargin)
    E = varargin{1};
    
    R = zeros(3,3,4);
    t = zeros(3,1,4);
    
    [U,S,V] = svd(E);
    B = [0 1 0;-1 0 0;0 0 1];
    
    sign1 = det(U*V');
    
    
    
    L = U*[0 -1 0;1 0 0;0 0 0]*U';
    M = -L;
    
    v1 = [L(3,2);L(1,3);L(2,1)];
    v2 = [M(3,2);M(1,3);M(2,1)];
    
    R(:,:,1) = sign1 * U * B * V';
    t(:,:,1) = v1./norm(v1);
    
    R(:,:,2) = sign1 * U * B' * V';
    t(:,:,2) = v2./norm(v2);
    
    R(:,:,3) = R(:,:,1);
    t(:,:,3) = t(:,:,2);
    
    R(:,:,4) = R(:,:,2);
    t(:,:,4) = t(:,:,1);
end