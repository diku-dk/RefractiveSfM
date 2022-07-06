
function ang = Rot2RPY(mat)

sinB = mat(3,1);
cosB = sqrt(mat(3,2)^2 + mat(3,3)^2);

if (abs(cosB) > 1e-15)
    ang(1) = atan2(mat(3,2), mat(3,3));
    ang(2) = atan2((-mat(3,1)), sqrt(mat(3,2)^2 + mat(3,3)^2));
    ang(3) = atan2(mat(2,1), mat(1,1));
else
    sinC = (mat(1,2) - mat(2,3)) / 2;
    cosC = (mat(2,2) + mat(1,3)) / 2;

    ang(1) = atan2(sinC, cosC); 

    ang(2) = M_PI/2; 
    ang(3) = 0;
    if (sinB < 0)
        for ii = 1:3
            ang(ii) = -ang(ii);
        end
    end
end

ang = ang * 180 / pi;


ang = ang(:);
