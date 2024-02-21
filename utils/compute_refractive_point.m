
function p = compute_refractive_point(rs,ds,n)
% rs is an array of rays;
% ds is an array of distance (layer);
% n is the common normal
% the following compute the refractive point.
% note: it assumes a common normal for refractive layers
    dots = dot(rs,repmat(n,1,size(rs,2)));% rs is 3xn, n is 3x1
    scales = ds ./ dots;% ds is 1xn
    p = rs .* scales;
end
