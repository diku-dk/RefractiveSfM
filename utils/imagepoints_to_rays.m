function ra = imagepoints_to_rays(K,q)
    if size(q,1) ~= 3
        % make it homogeneous
        q = [q;ones(1,size(q,2))];
    end
    % ray in general norm
    ra = inv(K) * q;
    % normalization to get the rays
    ra = ra./sqrt(ra(1,:).^2+ra(2,:).^2+ra(3,:).^2);
end