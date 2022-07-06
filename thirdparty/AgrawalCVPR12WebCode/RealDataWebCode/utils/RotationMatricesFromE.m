
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) MERL 2012
% CVPR 2012 Paper Title: A Theory of Multi-Layer Flat Refractive Geometry
% Author: Amit Agrawal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% get 4 rotation matrices from E matrix
% since we know axis direction, there are only 2 matrices. But there are
% two for -E also. Those 2 matrices need to be coupled with -s

function [Rsol,sSol] = RotationMatricesFromE(E,s,ne)

Rsol = zeros(3,3,4);
sSol = zeros(3,4);


% find rotation from E
[U,~,V] = svd(E);
%use Hartley matrices:
Wmat = [0 -1 0; 1 0 0; 0 0 1];
%Zmat = [0 1 0; -1 0 0; 0 0 0];
% get two rotations

R1 = U*Wmat*V';
R2 = U*Wmat'*V';

ncross = CrossMatrix(ne);

% find which solution is correct
e1 = ncross*R1 - E;
e1 = sum(sum(e1.^2));
e2 = ncross*R1 + E;
e2 = sum(sum(e2.^2));
if(e1 < e2)
    Rsol1 = R1;
else
    Rsol1 = -R1;
end

e3 = ncross*R2 - E;
e3 = sum(sum(e3.^2));
e4 = ncross*R2 + E;
e4 = sum(sum(e4.^2));
if(e3 < e4)
    Rsol2 = R2;
else
    Rsol2 = -R2;
end

Rsol(:,:,1) = Rsol1;
sSol(:,1) = s;

Rsol(:,:,2) = Rsol2;
sSol(:,2) = s;

Rsol(:,:,3) = -Rsol1;
sSol(:,3) = -s;

Rsol(:,:,4) = -Rsol2;
sSol(:,4) = -s;











