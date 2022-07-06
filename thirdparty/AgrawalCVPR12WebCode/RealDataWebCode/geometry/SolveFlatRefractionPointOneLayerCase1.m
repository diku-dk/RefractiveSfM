
% Solve Forward projection equation for Case 1 (single layer)
% d is the distance of medium from camera
% p is given 3D point
% n is the normal
% mu is refractive index (relative: so mu0/mu1)

% M is the 3D point on the layer closest to camera where the first
% refraction happens

function [M] = SolveFlatRefractionPointOneLayerCase1(d,n,mu,p)
    M = [0;0;1];

    %find POR the plane of refraction
    POR = cross(n,p);
    POR = POR/norm(POR);

    % [z1,z2] defines a coordinate system on POR
    % axis is away from the camera. z1 is along the axis
    z1 = -n;
    z1 = z1/norm(z1);

    % find z2
    z2 = cross(POR,z1);
    z2 = z2/norm(z2);

    % find the projection of given 3D point on POR
    v = p'*z1;
    u = p'*z2;

    %check
    %e = v*z1 + u*z2 - p

    %solve 4th degree equation on POR
    s1 = (mu - 1)*(mu + 1);
    s2 = (-2)*u*(mu - 1)*(mu + 1);
    s3 = d^2*mu^2 - d^2 + 2*d*v + mu^2*u^2 - u^2 - v^2;
    s4 = (-2)*d^2*mu^2*u;
    s5 = d^2*mu^2*u^2;

    sol = roots([s1;s2;s3;s4;s5]);

    % find real solutions
    idx = find(abs(imag(sol)) < 1e-6);
    if(isempty(idx))
        disp('no solution');
        return
    end

    sol1 = sol(idx);
    nn = size(sol1,1);

    % we need to find the correct solution out of all real solutions
    Normal = [0;-1];

    for ii = 1:nn

        x = sol1(ii,1);
        vi = [x;d];

        [v2,~,~] = RefractedRay(vi,Normal,1,mu);
        % this is actually the ray: point on POR - refractive point = refrative ray
        vrd = [u;v] - vi;

        e = abs(vrd(1)*v2(2) - vrd(2)*v2(1));

        if(e < 1e-4)
            M = x*z2 + d*z1;
            return
        end
    end
end




