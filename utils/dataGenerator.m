function varargout = dataGenerator(N,R1,n1,t1,refra_n1,refra_n2,varargin)
    if nargin == 6
        std = 0*1e-6;
    else
        std = varargin{1};
    end
    
    if nargin <= 7
        K = eye(3);
    else
        K = varargin{2};
    end
    
    d1 = abs(dot(n1,t1));  
    miu = refra_n2/refra_n1;%
    miu = 1/miu;
    for i = 1:N        
        while true
%         f1 = [rand(2,1)*0.3+0.1;1];x1 = f1./norm(f1);
        f1 = [rand(1)*1919+1;rand(1)*1079+1;1];x1 = inv(K)*f1;
        x1 = x1 ./ norm(x1);
        q1 = R1'*x1;q1 = q1 ./ norm(q1);
        ref_A = -d1 / dot(q1,n1) * q1 + t1;
        i1 = q1; % incident ray 1
        nn = -n1;
        tr1 = miu * i1 + sqrt(1-miu^2*(1-(dot(nn,i1))^2)) * nn - miu * dot(nn,i1) * nn;% transmitted ray 1
        tr1 = tr1./norm(tr1);
        if ~isreal(tr1)
            continue;
        end
        
        Q = ref_A + tr1 * (( 0.5 + 0.5 * rand(1) )/abs(tr1(3)));%+ 0.3 * rand(1) %+ 0.5 * rand(1) + 0.5 * rand(1)
        
        break
        end
        f1(1:2) = f1(1:2) + std.*randn(2,1);%q1 = q1 ./ norm(q1);
        q1s(:,i) = f1;
        Q1s(:,i) = Q;
        ref_As(:,i) = ref_A;
        plane_norm = cross(tr1,nn);plane_norm = plane_norm ./ norm(plane_norm);
        d = -dot(ref_A,plane_norm);
        planes(:,i) = [plane_norm;d];
        tr1s(:,i)=tr1;
    end
    varargout{1} = q1s(1:2,:);
%     figure
%     plot(q1s(1,:),q1s(2,:),'r.');
%     xlim([1,1920]);
%     ylim([1,1080]);
    varargout{2} = Q1s;
    varargout{3} = ref_As;
    varargout{4} = planes;
    varargout{5} = tr1s;
end