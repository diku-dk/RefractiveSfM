function [Rv,Cv,varargout] = find_virtual_camera_coordiante(ps,rws,n)
% this function computes the virtual camera coordiante system for RSFM as
% described in the paper:
% Refractive Structure-from-Motion on Underwater Images
    n = n./norm(n);
    Rv = cell(size(rws,2),1);
    Cv = cell(size(rws,2),1);
    tv = cell(size(rws,2),1);
    for i = 1:size(rws,2)
        rw = rws(:,i); p = ps(:,i);
        rw = rw./norm(rw);
        A = [n -rw];b = p;
        res = A\b;
        Cv{i} = p + res(2)*rw;
        rotaxis = cross([0;0;1],n);
        rotaxis = rotaxis./norm(rotaxis);
        rotangle = acos(dot([0;0;1],n));
        if abs(rotangle) < 1e-6
            Rv{i} = eye(3);
        else
            Rv{i} = expSO3(rotaxis*rotangle)';
        end
        tv{i} = -Rv{i}'*Cv{i};
    end
    varargout{1} = tv;
end