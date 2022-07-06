function [Cost] = CostFuncOneLayersCase1(R,C,n,d,rayDnoise,XYZ,KK,mu)
    N = size(XYZ,2);
    Cost = zeros(2*N,1);
    ps = R'*(XYZ - C);
    for ii = 1:N
        p = ps(:,ii);
        M = SolveFlatRefractionPointOneLayerCase1(d,n,mu,p);
        M = KK*M;
        px1 = M(1)/M(3);
        py1 = M(2)/M(3);
        px = rayDnoise(1,ii);
        py = rayDnoise(2,ii);
        Cost(2*ii-1,1) = px-px1;
        Cost(2*ii,1) = py-py1;
    end
end