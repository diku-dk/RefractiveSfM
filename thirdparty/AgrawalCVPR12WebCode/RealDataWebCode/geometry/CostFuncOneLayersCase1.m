
% Find vector of reprojection errors for non-linear refinement
% for one layer Case 1

function [C] = CostFuncOneLayersCase1(x,rayDnoise,XYZ,KK,mu)
    x = x(:);

    n = x(1:2);
    w = x(3:5);
    t = x(6:8);
    tau = x(9);
    n = [n ; -1];



    n = n/norm(n);

    R = RPY2Rot(w);

    N = size(XYZ,2);

    C = zeros(2*N,1);

    for ii = 1:N
        p = R*XYZ(:,ii) + t;


        M = SolveFlatRefractionPointOneLayerCase1(tau,n,mu,p);
        M = KK*M;
        px1 = M(1)/M(3);
        py1 = M(2)/M(3);
        px = rayDnoise(1,ii);
        py = rayDnoise(2,ii);
        C(2*ii-1,1) = px-px1;
        C(2*ii,1) = py-py1;
    end
end









