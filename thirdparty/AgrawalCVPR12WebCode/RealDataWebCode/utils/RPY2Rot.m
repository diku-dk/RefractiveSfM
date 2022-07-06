


%function mat = RPY2Rot(roll, pitch, yaw)

function mat = RPY2Rot(x)

roll = x(1);
pitch = x(2);
yaw = x(3);

r = roll * pi / 180;
p = pitch * pi / 180;
y = yaw * pi / 180;

rollMat = [ 1 0 0 0;
            0 cos(r) -sin(r) 0;
            0 sin(r) cos(r) 0;
            0 0 0 1];
        
pitchMat = [cos(p) 0 sin(p) 0;
            0 1 0 0;
            -sin(p) 0 cos(p) 0;
            0 0 0 1];
            
yawMat = [  cos(y) -sin(y) 0 0;
            sin(y) cos(y) 0 0;
            0 0 1 0;
            0 0 0 1];
        
mat = yawMat * pitchMat * rollMat;


mat = mat(1:3,1:3);

