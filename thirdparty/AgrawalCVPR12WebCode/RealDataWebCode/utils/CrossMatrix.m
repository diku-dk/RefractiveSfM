
function [ncross] = CrossMatrix(n)

ncross = [0 -n(3) n(2) ; n(3) 0 -n(1) ; -n(2) n(1) 0];
