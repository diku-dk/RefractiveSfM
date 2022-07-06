


function w = R2w(R)


theta = acos((trace(R)-1)/2);

if(abs(theta*180/pi) < 1e-6)
    % rotation angle is 0
    w = zeros(3,1);
    return
end


if(abs(theta-pi) < 1e-3)
    % angle of rotation is 180 degrees
    w = [pi ; 0; 0 ];
    return
end

s = sin(theta);

w = [ R(3,2) - R(2,3) ; R(1,3) - R(3,1) ; R(2,1) - R(1,2)] / s;

if(norm(w) > 1e-6)
    w = theta*w/norm(w);
end





% % find the eigen-vector which has eigenvalue 1
% [v,d] = eig(R);
%
% d = diag(d);
%
% e = abs(d-1);
% [~,id] = min(e);
%
% w = v(:,id);
%
% % trace = 1+2cos(theta)
% mag = abs(acos((trace(R) - 1)/2));
%
% w = mag*w/norm(w);
%
%






