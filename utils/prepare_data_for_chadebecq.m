function [NrP, Chats] = prepare_data_for_chadebecq(lambda,d1c,d2c,NN,q2s)
    Ss3=[1 0 0 0 0 0 0 0 0;...
        0 1 0 1 0 0 0 0 0;...
        0 0 0 0 1 0 0 0 0;...
        0 0 1 0 0 0 1 0 0;...
        0 0 0 0 0 1 0 1 0;...
        0 0 0 0 0 0 0 0 1];
    Ds3=(diag(sum(Ss3,2)));
    invDs3 = inv(Ds3);
    lambda2 = lambda^2;
    fLambda = @(t) (invDs3*Ss3*kron(t,t)*Ss3');
    tb = [0 0 1;0 d1c 0;-d1c 0 0]';
    ta = [1 0 0;0 1 0;0 0 0]';
    tmp1 = [(1-lambda2)*fLambda(tb) zeros(6,6);...
            (lambda2)*fLambda(tb) -(lambda2)*fLambda(ta)]';
    Ds3bar = blkdiag(Ds3,Ds3);

    fhat6 = @(q) ([q(1)^2 q(1)*q(2) q(2)^2 q(1)*q(3) q(2)*q(3) q(3)^2 q(1)*q(4) q(2)*q(4) q(3)*q(4) q(4)^2 q(1)*q(5) q(2)*q(5) q(3)*q(5) q(4)*q(5) q(5)^2 q(1)*q(6) q(2)*q(6) q(3)*q(6) q(4)*q(6) q(5)*q(6) q(6)^2]');
    fhat = @(q) ([q(1)*q(1);q(1)*q(2);q(2)*q(2);q(1)*q(3);q(2)*q(3);q(3)*q(3)]);

    rP = tmp1*Ds3bar;

    tu = [0 0 0;1 0 0;0 1 0];
    tv = [0 d2c 0;-d2c 0 0;0 0 1];
    Ss = zeros(21,36);
    ind = [1,1;2,2;2,7;3,8;4 3;4 13;5 9;5 14;6 15;7 4;7 19;8 10;8 20;9 16;9 21;10 22;11 5;11 25;12 11;12 26;13 17;13 27;14 23;14 28;15 29;16 6;16 31;17 12;17 32;18 18;18 33;19 24;19 34;20 30;20 35;21 36];
    Ss(ind(:,1)+(ind(:,2)-1)*21) = 1;
    Ds=(diag(sum(Ss,2)));

    indN = [1 1; 2 2;3 3;4 4;5 5;6 6;7 10;8 14;9 15;10 19;11 20;12 21];
    N = zeros(12,21);
    N(indN(:,1) + (indN(:,2)-1)*12) = 1;
    
    % firstly check peter strum 
    Chats = cell(NN,1);
    for i = 1:NN
        gemma = sqrt(1-lambda^2+lambda^2*q2s(3,i)^2);
        C = [lambda*tu zeros(3,3);zeros(3,3) gemma*tv];
        Chat = inv(Ds)*Ss*kron(C,C)*Ss';
        Chats{i} = Chat;
    end
    NrP = N'*rP;
end