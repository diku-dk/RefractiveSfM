


% Find the subspace using 8 points and give it to 5-point algorithm to
% estimate the coeffs

% Directly taken from 5 point implementation by Nister/Stewanius

% done

function [SOLS] = calibrated_fivepoint_Amit(EE)


A = calibrated_fivepoint_helper(EE) ;

KK = A(:,1:10);

if(rank(KK) < 10)
    SOLS = Inf*ones(1,4);
    return
end


A = A(:,1:10)\A(:,11:20);
M = -A([1 2 3 5 6 8], :);

M(7,1) = 1;
M(8,2) = 1;
M(9,4) = 1;
M(10,7) = 1;

[V,~] = eig(M );
SOLS =   V(7:9,:)./(ones(3,1)*V(10,:));

SOLS = [SOLS ; ones(1,10 ) ];

idx = not(imag(SOLS(1,:)));
SOLS = SOLS(:,idx);

SOLS = SOLS';

