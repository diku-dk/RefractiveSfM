


%R = I\cos\theta + \sin\theta[\mathbf u]_{\times} + (1-\cos\theta)\mathbf{u}\otimes\mathbf{u},


function R = w2R(w)


theta = norm(w);

if(abs(theta*180/pi) < 1e-6)
    R = eye(3);
    return
end

    
w = w/theta;

c = cos(theta);
s = sin(theta);

wcross = [0 -w(3) w(2) ; w(3) 0 -w(1) ; -w(2) w(1) 0];

R = eye(3)*c + s*wcross + (1-c)*(w*w');


