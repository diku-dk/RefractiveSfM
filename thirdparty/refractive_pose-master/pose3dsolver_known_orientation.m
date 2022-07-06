function s = translation3dsolver(Rcam,X,x,N,W,nrat)


% Transform coordinates so that the refractive plane coincides with z = 0.
[U,~,~] = svd(N);
R = [U(:,2)*det(U) U(:,[3 1])];
X = R'*(X-W(:)*[1 1]);
u = [x(1:2,:); 1 1];
u = bsxfun(@rdivide,u,sqrt(sum(u.^2,1)));
u = R'*Rcam'*u;

% Compute equation coefficients
C = eqcoeff(u(:,1),u(:,2),X(:,1),X(:,2));

if nargin<6 % Unknown refractive index ratio

	coef = [C(7)*C(4) - C(1)*C(10),...
			C(7)*C(5) - C(1)*C(11) - C(2)*C(10) + C(8)*C(4),...
			C(7)*C(6) - C(1)*C(12) - C(2)*C(11) + C(8)*C(5) - C(3)*C(10) + C(9)*C(4),...
			C(8)*C(6) - C(2)*C(12) - C(3)*C(11) + C(9)*C(5),...
			C(9)*C(6) - C(3)*C(12)];
		
	tz = roots(coef);

	r2 = -(C(4)*tz.^2 + C(5)*tz + C(6))./(C(1)*tz.^2 + C(2)*tz + C(3));

	txy = [u(2,1) -u(1,1); u(2,2) -u(1,2)]\[u(2,1)*X(1,1) - u(1,1)*X(2,1); u(2,2)*X(1,2) - u(1,2)*X(2,2)];
	
	% Transform to original coordinates
	t = [repmat(txy,1,4); tz(:).'];
	t = R*t + W(:)*[1 1 1 1];
	s = [t; sqrt(r2(:).')];
	
else
	% Here we can choose which quadratic equation to use. With perfect data
	% these will have one common root corresponding to the true translation
	% while giving different non-physical alternate solutions.
	coef = [C(1)*nrat^2 + C(4), C(2)*nrat^2 + C(5), C(3)*nrat^2 + C(6)];
% 	coef = [C(7)*nrat^2 + C(10), C(8)*nrat^2 + C(11), C(9)*nrat^2 + C(12)];


	tz = roots(coef);
	
	txy = [u(2,1) -u(1,1); u(2,2) -u(1,2)]\[u(2,1)*X(1,1) - u(1,1)*X(2,1); u(2,2)*X(1,2) - u(1,2)*X(2,2)];
	
	% Transform to original coordinates
	t = [repmat(txy,1,2); tz(:).'];
	s = R*t + W(:)*[1 1];
	
end




function A0 = eqcoeff(u1,u2,X1,X2)

u11 = u1(1); u21 = u1(2); u31 = u1(3);
u12 = u2(1); u22 = u2(2); u32 = u2(3);
X11 = X1(1); X21 = X1(2); X31 = X1(3);
X12 = X2(1); X22 = X2(2); X32 = X2(3);

A0 = zeros(1,12);

t2 = u21*u21;
t3 = u11*u22;
t7 = u12*u21;
t4 = t3-t7;
t5 = 1.0/t4;
t6 = X11*u31;
t8 = X11*u12*u21;
t9 = X22*u11*u12;
t24 = X21*u11*u12;
t25 = X12*u11*u22;
t10 = t8+t9-t24-t25;
t11 = t5*t10*u31;
t12 = t6+t11;
t13 = X21*u31;
t14 = X11*u21*u22;
t15 = X22*u12*u21;
t19 = X21*u11*u22;
t20 = X12*u21*u22;
t16 = t14+t15-t19-t20;
t17 = t5*t16*u31;
t18 = t13+t17;
t21 = t18*u21*2.0;
t22 = t18*t18;
t23 = u22*u22;
t26 = X12*u32;
t27 = t5*t10*u32;
t28 = t26+t27;
t29 = X22*u32;
t30 = t5*t16*u32;
t31 = t29+t30;
t32 = t31*u22*2.0;
t33 = t31*t31;
A0(1,1) = t2*(t2+u11*u11);
A0(1,2) = t2*(t21+t12*u11*2.0);
A0(1,3) = t2*(t22+t12*t12+(X31*X31)*(u31*u31));
A0(1,4) = -t2;
A0(1,5) = -t21;
A0(1,6) = -t22;
A0(1,7) = t23*(t23+u12*u12);
A0(1,8) = t23*(t32+t28*u12*2.0);
A0(1,9) = t23*(t33+t28*t28+(X32*X32)*(u32*u32));
A0(1,10) = -t23;
A0(1,11) = -t32;
A0(1,12) = -t33;

