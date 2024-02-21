
function rg = compute_refractive_ray(na,nb,ra,n)
% na is where the input ray exists
% nb is where the output ray exists
% ra is the input ray
% n is the normal
% the following equation works well in the camera coordinate system since
% the origin is 0 0 0
    miu = na/nb;
    nn = repmat(n,1,size(ra,2));
    rg = miu * ra + (-miu * dot(nn,ra) + sqrt(1-miu^2*(1-(dot(nn,ra).^2)))).*n;
    rg = rg./sqrt(rg(1,:).^2+rg(2,:).^2+rg(3,:).^2);
end