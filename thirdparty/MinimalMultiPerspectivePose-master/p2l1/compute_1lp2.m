function [sR,sT,sdepth1,sdepth2,aR1,aR2,at2] = compute_1lp2(adapter)
% AUTORIGHTS
% ---------------------------------------------------------
% Copyright (c) 2018, Pedro Miraldo
% 
% This file is part of the Minimal Multi-Perspective Pose
% code and is available under the terms of the MIT License
% provided in LICENSE. Please retain this notice and
% LICENSE if you use this file (or any portion of it)
% in your project.
% ---------------------------------------------------------

aR1 = 0;
aR2 = 0;
at2 = 0;

% get data from the adapter
% camera info
C1 = adapter.camera.c1;
C2 = adapter.camera.c2;
D1 = adapter.camera.d1;
D2 = adapter.camera.d2;
N1 = adapter.camera.l1.normal;
% world info
P1  = adapter.world.p1;
P2  = adapter.world.p2;
LP1 = adapter.world.l1.p1;
LP2 = adapter.world.l1.p2;

% get the right coordinate systems
N11 = N1(1); N12 = N1(2); N13 = N1(3);

% Get the correct camera coordinates frame (a Rotation is enough)
e12_ =  N1(3); e13_ = -N1(2);
e11__ = N1(3); e13__ = -N1(1);

ne1_ = e12_*e12_ + e13_*e13_;
ne1__ = e11__*e11__ + e13__*e13__;
if ne1_ >= ne1__ ;
    e11 = 0;
    e12 = e12_/(ne1_)^(1/2);
    e13 = e13_/(ne1_)^(1/2);
else
    e11 = e11__/(ne1__)^(1/2);
    e12 = 0;
    e13 = e13__/(ne1__)^(1/2);
end

e21 = -N13*e12+N12*e13;
e22 = +N13*e11-N11*e13;
e23 = -N12*e11+N11*e12;

e31 = N11;
e32 = N12;
e33 = N13;

% Apply the transformation to the data
% D1
d11 = e11*D1(1) + e12*D1(2) + e13*D1(3);
d12 = e21*D1(1) + e22*D1(2) + e23*D1(3);
d13 = e31*D1(1) + e32*D1(2) + e33*D1(3);
% D2
d21 = e11*D2(1) + e12*D2(2) + e13*D2(3);
d22 = e21*D2(1) + e22*D2(2) + e23*D2(3);
d23 = e31*D2(1) + e32*D2(2) + e33*D2(3);
% C1
c11 = e11*C1(1) + e12*C1(2) + e13*C1(3);
c12 = e21*C1(1) + e22*C1(2) + e23*C1(3);
c13 = e31*C1(1) + e32*C1(2) + e33*C1(3);
% C2
c21 = e11*C2(1) + e12*C2(2) + e13*C2(3);
c22 = e21*C2(1) + e22*C2(2) + e23*C2(3);
c23 = e31*C2(1) + e32*C2(2) + e33*C2(3);

% Transformation to the world data
DL1 = LP2-LP1;
f21 = DL1(1);
f22 = DL1(2);
f23 = DL1(3);
nf2 = (f21*f21+f22*f22+f23*f23)^(1/2);
f21 = f21/nf2; f22 = f22/nf2; f23 = f23/nf2;

n1_ = P1(1)-LP1(1);
n2_ = P1(2)-LP1(2);
n3_ = P1(3)-LP1(3);
f11 = -f23*n2_+f22*n3_;
f12 = +f23*n1_-f21*n3_;
f13 = -f22*n1_+f21*n2_;
nf1 = (f11*f11+f12*f12+f13*f13)^(1/2);
f11 = f11/nf1;
f12 = f12/nf1;
f13 = f13/nf1;

f31 = -f13*f22+f12*f23;
f32 = +f13*f21-f11*f23;
f33 = -f12*f21+f11*f22;


p12_ = f21*P1(1) + f22*P1(2) + f23*P1(3);
LP12_ = f21*LP1(1) + f22*LP1(2) + f23*LP1(3);
l = -LP12_ + p12_;

t1_ = -(f11*LP1(1) + f12*LP1(2) + f13*LP1(3));
t2_ = -(f21*LP1(1) + f22*LP1(2) + f23*LP1(3) + l);
t3_ = -(f31*LP1(1) + f32*LP1(2) + f33*LP1(3));

p13 = f31*P1(1) + f32*P1(2) + f33*P1(3) + t3_;

p21 = f11*P2(1) + f12*P2(2) + f13*P2(3) + t1_;
p22 = f21*P2(1) + f22*P2(2) + f23*P2(3) + t2_;
p23 = f31*P2(1) + f32*P2(2) + f33*P2(3) + t3_;

% compute the pose
cc11 = -d13^2*(p21^2 + p23^2);
cc12 = 2*d13*d23*p13*p23;
cc13 = -2*d13*(c13*p21^2 + c13*p23^2 - c23*p13*p23);
cc14 = -d23^2*p13^2;
cc15 = 2*d23*p13*(c13*p23 - c23*p13);
cc16 = - c13^2*p21^2 - c13^2*p23^2 + 2*c13*c23*p13*p23 - c23^2*p13^2 + p13^2*p21^2;

cc21 = - d11^2*p13^2*p21^2 - d12^2*p13^2*p21^2 + d13^2*p13^2*p23^2 - 2*d13^2*p13*p21^2*p23 - 2*d13^2*p13*p23^3 + d13^2*p21^4 + 2*d13^2*p21^2*p23^2 + d13^2*p23^4;
cc22 = 2*d11*d21*p13^2*p21^2 - 2*d13*d23*p13^3*p23 - 2*d13*d23*p13*p23^3 + 2*d12*d22*p13^2*p21^2 + 2*d13*d23*p13^2*p21^2 + 4*d13*d23*p13^2*p23^2 - 2*d13*d23*p13*p21^2*p23;
cc23 = 2*c13*d13*p21^4 + 2*c13*d13*p23^4 - 4*c13*d13*p13*p23^3 - 2*c23*d13*p13*p23^3 - 2*c23*d13*p13^3*p23 - 2*c11*d11*p13^2*p21^2 - 2*c12*d12*p13^2*p21^2 + 2*c13*d13*p13^2*p23^2 + 2*c21*d11*p13^2*p21^2 + 2*c22*d12*p13^2*p21^2 + 4*c13*d13*p21^2*p23^2 + 2*c23*d13*p13^2*p21^2 + 4*c23*d13*p13^2*p23^2 - 4*c13*d13*p13*p21^2*p23 - 2*c23*d13*p13*p21^2*p23;
cc24 = -p13^2*(d21^2*p21^2 + d22^2*p21^2 - d23^2*p13^2 + 2*d23^2*p13*p23 - d23^2*p23^2);
cc25 = 2*c23*d23*p13^4 - 2*c13*d23*p13*p23^3 - 2*c13*d23*p13^3*p23 - 4*c23*d23*p13^3*p23 + 2*c11*d21*p13^2*p21^2 + 2*c12*d22*p13^2*p21^2 + 2*c13*d23*p13^2*p21^2 + 4*c13*d23*p13^2*p23^2 - 2*c21*d21*p13^2*p21^2 - 2*c22*d22*p13^2*p21^2 + 2*c23*d23*p13^2*p23^2 - 2*c13*d23*p13*p21^2*p23;
cc26 = - c11^2*p13^2*p21^2 + 2*c11*c21*p13^2*p21^2 - c12^2*p13^2*p21^2 + 2*c12*c22*p13^2*p21^2 + c13^2*p13^2*p23^2 - 2*c13^2*p13*p21^2*p23 - 2*c13^2*p13*p23^3 + c13^2*p21^4 + 2*c13^2*p21^2*p23^2 + c13^2*p23^4 - 2*c13*c23*p13^3*p23 + 2*c13*c23*p13^2*p21^2 + 4*c13*c23*p13^2*p23^2 - 2*c13*c23*p13*p21^2*p23 - 2*c13*c23*p13*p23^3 - c21^2*p13^2*p21^2 - c22^2*p13^2*p21^2 + c23^2*p13^4 - 2*c23^2*p13^3*p23 + c23^2*p13^2*p23^2 + p13^2*p21^2*p22^2;

fc1 = cc11^2 + (cc12^2*cc21)/cc24 + (cc14^2*cc21^2)/cc24^2 - (cc11*cc12*cc22)/cc24 - (2*cc11*cc14*cc21)/cc24 + (cc11*cc14*cc22^2)/cc24^2 - (cc12*cc14*cc21*cc22)/cc24^2;
fc2 = 2*cc11*cc13 + (cc12^2*cc23)/cc24 - (cc12*cc13*cc22)/cc24 - (cc11*cc12*cc25)/cc24 - (2*cc11*cc14*cc23)/cc24 - (cc11*cc15*cc22)/cc24 + (2*cc12*cc15*cc21)/cc24 - (2*cc13*cc14*cc21)/cc24 + (cc13*cc14*cc22^2)/cc24^2 + (2*cc14^2*cc21*cc23)/cc24^2 - (cc12*cc14*cc22*cc23)/cc24^2 + (2*cc11*cc14*cc22*cc25)/cc24^2 - (cc12*cc14*cc21*cc25)/cc24^2 - (cc14*cc15*cc21*cc22)/cc24^2;
fc3 = 2*cc11*cc16 + cc13^2 + (cc15^2*cc21)/cc24 + (cc12^2*cc26)/cc24 + (cc14^2*cc23^2)/cc24^2 - (cc12*cc13*cc25)/cc24 + (2*cc12*cc15*cc23)/cc24 - (cc12*cc16*cc22)/cc24 - (2*cc13*cc14*cc23)/cc24 - (cc13*cc15*cc22)/cc24 - (2*cc11*cc14*cc26)/cc24 - (cc11*cc15*cc25)/cc24 - (2*cc14*cc16*cc21)/cc24 + (cc11*cc14*cc25^2)/cc24^2 + (cc14*cc16*cc22^2)/cc24^2 + (2*cc14^2*cc21*cc26)/cc24^2 - (cc12*cc14*cc22*cc26)/cc24^2 - (cc12*cc14*cc23*cc25)/cc24^2 + (2*cc13*cc14*cc22*cc25)/cc24^2 - (cc14*cc15*cc22*cc23)/cc24^2 - (cc14*cc15*cc21*cc25)/cc24^2;
fc4 = 2*cc13*cc16 + (cc15^2*cc23)/cc24 + (2*cc12*cc15*cc26)/cc24 - (cc12*cc16*cc25)/cc24 - (2*cc13*cc14*cc26)/cc24 - (cc13*cc15*cc25)/cc24 - (2*cc14*cc16*cc23)/cc24 - (cc15*cc16*cc22)/cc24 + (cc13*cc14*cc25^2)/cc24^2 + (2*cc14^2*cc23*cc26)/cc24^2 - (cc12*cc14*cc25*cc26)/cc24^2 - (cc14*cc15*cc22*cc26)/cc24^2 - (cc14*cc15*cc23*cc25)/cc24^2 + (2*cc14*cc16*cc22*cc25)/cc24^2;
fc5 = cc16^2 + (cc15^2*cc26)/cc24 + (cc14^2*cc26^2)/cc24^2 - (2*cc14*cc16*cc26)/cc24 - (cc15*cc16*cc25)/cc24 + (cc14*cc16*cc25^2)/cc24^2 - (cc14*cc15*cc25*cc26)/cc24^2;

%% ferrari formula

b = fc2/fc1; c = fc3/fc1; d = fc4/fc1;
e = fc5/fc1; a = fc1/fc1;

f = c - ((3*b*b)/8);
g = d + ((b*b*b)/8) - ((b*c)/2);
h = e - ((3*b*b*b*b)/256) + (((b*b)*c)/16) -((b*d)/4);

%Cubic polynomial
a2 = 1;
b2 = f/2;
c2 = ((f*f - 4*h)/16);
d2 = -(g*g)/64;

f2 = (((3*c2)/a2) - ((b2*b2)/(a2*a2)))/3;
g2 = (((2*(b2*b2*b2))/(a2*a2*a2)) - ((9*b2*c2)/(a2*a2)) + ((27*d2)/a2))/27;
h2 = ((g2*g2)/4) + ((f2*f2*f2)/27);

i2 = sqrt(((g2*g2)/4) - h2);
j2 = i2^(1/3);
k2 = acos(-(g2/(2*i2)));

l2 = -j2;
m2 = cos(k2/3);
n2 = sqrt(3)*sin(k2/3);
p2 = -(b2/(3*a2));

x1 = 2*j2 * cos(k2/3) - (b2/(3*a2));
x2 = l2 * (m2 + n2) + p2;
x3 = l2 * (m2 - n2) + p2;

if(x1~= 0 && x2~=0)
    p = sqrt(x1);
    q = sqrt(x2);
elseif (x1 ~= 0 && x3 ~= 0)
    p = sqrt(x1);
    q = sqrt(x3);
elseif (x2 ~= 0 && x3 ~=0)
    p = sqrt(x2);
    q = sqrt(x3);
else
    disp('Error');
    return;
end

r = -g/(8*p*q);
s = b/(4*a);

a1 = p + q + r -s;
a2 = p - q - r - s;
a3 = -p + q - r -s;
a4 = -p - q + r -s;

c = [a1,a2,a3,a4];

% c = roots([fc1,fc2,fc3,fc4,fc5]);

sdepth1 = real(c(abs(imag(c))<1e-4)');
if isempty(sdepth1)
    sR = [];sT= [];
    sdepth1 = []; sdepth2 = [];
    return;
end

% set the output
sdepth2 = zeros(1,numel(sdepth1));
sR = zeros(3,3*numel(sdepth1));
sT = zeros(3,numel(sdepth1));

for i = 1 : numel(sdepth1)
    
    alpha1 = sdepth1(i);
    
    vp = (alpha1^2*cc22^2 - 4*cc21*cc24*alpha1^2 + 2*alpha1*cc22*cc25 - 4*cc23*cc24*alpha1 + cc25^2 - 4*cc24*cc26);
    p2 = -(cc25 + alpha1*cc22);
    p3 = (2*cc24);
    
    l1 = alpha1;
    l2_ = (p2 + vp^(1/2))/p3;
    l2__ = (p2 - vp^(1/2))/p3;
    
    % New hyp
    if abs(imag(l2_))<10e-4
        l2 = l2_;
        
        % get R hyp
        frac = c13^2*p21^2 + c13^2*p23^2 - 2*c13*c23*p13*p23 + 2*c13*d13*l1*p21^2 + 2*c13*d13*l1*p23^2 - 2*c13*d23*l2*p13*p23 + c23^2*p13^2 - 2*c23*d13*l1*p13*p23 + 2*c23*d23*l2*p13^2 + d13^2*l1^2*p21^2 + d13^2*l1^2*p23^2 - 2*d13*d23*l1*l2*p13*p23 + d23^2*l2^2*p13^2;
        c1 = (c13*p13*p21^2 + d13*l1*p13*p21^2)/frac;
        s1 = (p13*p21*(c13*p23 - c23*p13 + d13*l1*p23 - d23*l2*p13))/frac;
        frac = p13*p21*(c11^2 - 2*c11*c21 + 2*c11*d11*l1 - 2*c11*d21*l2 + c12^2 - 2*c12*c22 + 2*c12*d12*l1 - 2*c12*d22*l2 + c21^2 - 2*c21*d11*l1 + 2*c21*d21*l2 + c22^2 - 2*c22*d12*l1 + 2*c22*d22*l2 + d11^2*l1^2 - 2*d11*d21*l1*l2 + d12^2*l1^2 - 2*d12*d22*l1*l2 + d21^2*l2^2 + d22^2*l2^2);
        c2 = -(c11*c13*p21^2 + c11*c13*p23^2 + c11*c23*p13^2 - c13*c21*p21^2 - c13*c21*p23^2 - c21*c23*p13^2 - c11*c13*p13*p23 - c11*c23*p13*p23 + c13*c21*p13*p23 + c21*c23*p13*p23 + c12*p13*p21*p22 - c22*p13*p21*p22 + c11*d13*l1*p21^2 + c13*d11*l1*p21^2 + c11*d13*l1*p23^2 + c13*d11*l1*p23^2 + c23*d11*l1*p13^2 + c11*d23*l2*p13^2 - c21*d13*l1*p21^2 - c13*d21*l2*p21^2 - c21*d13*l1*p23^2 - c13*d21*l2*p23^2 - c21*d23*l2*p13^2 - c23*d21*l2*p13^2 + d11*d13*l1^2*p21^2 + d11*d13*l1^2*p23^2 - d21*d23*l2^2*p13^2 - c11*d13*l1*p13*p23 - c13*d11*l1*p13*p23 + c21*d13*l1*p13*p23 - c23*d11*l1*p13*p23 - c11*d23*l2*p13*p23 + c13*d21*l2*p13*p23 + c21*d23*l2*p13*p23 + c23*d21*l2*p13*p23 + d12*l1*p13*p21*p22 - d22*l2*p13*p21*p22 + d11*d23*l1*l2*p13^2 - d13*d21*l1*l2*p21^2 - d13*d21*l1*l2*p23^2 - d11*d13*l1^2*p13*p23 + d21*d23*l2^2*p13*p23 - d11*d23*l1*l2*p13*p23 + d13*d21*l1*l2*p13*p23)/frac;
        s2 = -(c12*c13*p21^2 + c12*c13*p23^2 + c12*c23*p13^2 - c13*c22*p21^2 - c13*c22*p23^2 - c22*c23*p13^2 - c12*c13*p13*p23 - c12*c23*p13*p23 + c13*c22*p13*p23 + c22*c23*p13*p23 - c11*p13*p21*p22 + c21*p13*p21*p22 + c12*d13*l1*p21^2 + c13*d12*l1*p21^2 + c12*d13*l1*p23^2 + c13*d12*l1*p23^2 + c23*d12*l1*p13^2 + c12*d23*l2*p13^2 - c22*d13*l1*p21^2 - c13*d22*l2*p21^2 - c22*d13*l1*p23^2 - c13*d22*l2*p23^2 - c22*d23*l2*p13^2 - c23*d22*l2*p13^2 + d12*d13*l1^2*p21^2 + d12*d13*l1^2*p23^2 - d22*d23*l2^2*p13^2 - c12*d13*l1*p13*p23 - c13*d12*l1*p13*p23 + c22*d13*l1*p13*p23 - c23*d12*l1*p13*p23 - c12*d23*l2*p13*p23 + c13*d22*l2*p13*p23 + c22*d23*l2*p13*p23 + c23*d22*l2*p13*p23 - d11*l1*p13*p21*p22 + d21*l2*p13*p21*p22 + d12*d23*l1*l2*p13^2 - d13*d22*l1*l2*p21^2 - d13*d22*l1*l2*p23^2 - d12*d13*l1^2*p13*p23 + d22*d23*l2^2*p13*p23 - d12*d23*l1*l2*p13*p23 + d13*d22*l1*l2*p13*p23)/frac;
        
        % get T hyp
        t1 = s1*(c13 + d13*l1) - c1*c2*(c11 + d11*l1) - c1*s2*(c12 + d12*l1);
        t2 = s2*(c11 + d11*l1) - c2*(c12 + d12*l1);
        t3 = p13 - c1*(c13 + d13*l1) - c2*s1*(c11 + d11*l1) - s1*s2*(c12 + d12*l1);
        
        if abs(c1*c1 + s1*s1 - 1) < 10e-2 && abs(c2*c2 + s2*s2 - 1) < 10e-2
            sdepth1(i) = l1;
            sdepth2(i) = l2;
            sR(:,(i-1)*3+1:i*3) = [ e11*(c1*c2*f11 - f21*s2 + c2*f31*s1) + e31*(c1*f31 - f11*s1) + e21*(c2*f21 + c1*f11*s2 + f31*s1*s2), e12*(c1*c2*f11 - f21*s2 + c2*f31*s1) + e32*(c1*f31 - f11*s1) + e22*(c2*f21 + c1*f11*s2 + f31*s1*s2), e13*(c1*c2*f11 - f21*s2 + c2*f31*s1) + e33*(c1*f31 - f11*s1) + e23*(c2*f21 + c1*f11*s2 + f31*s1*s2);
                e11*(c1*c2*f12 - f22*s2 + c2*f32*s1) + e31*(c1*f32 - f12*s1) + e21*(c2*f22 + c1*f12*s2 + f32*s1*s2), e12*(c1*c2*f12 - f22*s2 + c2*f32*s1) + e32*(c1*f32 - f12*s1) + e22*(c2*f22 + c1*f12*s2 + f32*s1*s2), e13*(c1*c2*f12 - f22*s2 + c2*f32*s1) + e33*(c1*f32 - f12*s1) + e23*(c2*f22 + c1*f12*s2 + f32*s1*s2);
                e11*(c1*c2*f13 - f23*s2 + c2*f33*s1) + e31*(c1*f33 - f13*s1) + e21*(c2*f23 + c1*f13*s2 + f33*s1*s2), e12*(c1*c2*f13 - f23*s2 + c2*f33*s1) + e32*(c1*f33 - f13*s1) + e22*(c2*f23 + c1*f13*s2 + f33*s1*s2), e13*(c1*c2*f13 - f23*s2 + c2*f33*s1) + e33*(c1*f33 - f13*s1) + e23*(c2*f23 + c1*f13*s2 + f33*s1*s2)];
            sT(:,i) =  [f11*t1 + f21*t2 + f31*t3 - f11*t1_ - f21*t2_ - f31*t3_
                f12*t1 + f22*t2 + f32*t3 - f12*t1_ - f22*t2_ - f32*t3_
                f13*t1 + f23*t2 + f33*t3 - f13*t1_ - f23*t2_ - f33*t3_];
        end
        
    end
    
    % New hyp
    if abs(imag(l2__))<10e-4
        
        l2 = l2__;
        
        % get R hyp
        frac = c13^2*p21^2 + c13^2*p23^2 - 2*c13*c23*p13*p23 + 2*c13*d13*l1*p21^2 + 2*c13*d13*l1*p23^2 - 2*c13*d23*l2*p13*p23 + c23^2*p13^2 - 2*c23*d13*l1*p13*p23 + 2*c23*d23*l2*p13^2 + d13^2*l1^2*p21^2 + d13^2*l1^2*p23^2 - 2*d13*d23*l1*l2*p13*p23 + d23^2*l2^2*p13^2;
        c1 = (c13*p13*p21^2 + d13*l1*p13*p21^2)/frac;
        s1 = (p13*p21*(c13*p23 - c23*p13 + d13*l1*p23 - d23*l2*p13))/frac;
        frac = p13*p21*(c11^2 - 2*c11*c21 + 2*c11*d11*l1 - 2*c11*d21*l2 + c12^2 - 2*c12*c22 + 2*c12*d12*l1 - 2*c12*d22*l2 + c21^2 - 2*c21*d11*l1 + 2*c21*d21*l2 + c22^2 - 2*c22*d12*l1 + 2*c22*d22*l2 + d11^2*l1^2 - 2*d11*d21*l1*l2 + d12^2*l1^2 - 2*d12*d22*l1*l2 + d21^2*l2^2 + d22^2*l2^2);
        c2 = -(c11*c13*p21^2 + c11*c13*p23^2 + c11*c23*p13^2 - c13*c21*p21^2 - c13*c21*p23^2 - c21*c23*p13^2 - c11*c13*p13*p23 - c11*c23*p13*p23 + c13*c21*p13*p23 + c21*c23*p13*p23 + c12*p13*p21*p22 - c22*p13*p21*p22 + c11*d13*l1*p21^2 + c13*d11*l1*p21^2 + c11*d13*l1*p23^2 + c13*d11*l1*p23^2 + c23*d11*l1*p13^2 + c11*d23*l2*p13^2 - c21*d13*l1*p21^2 - c13*d21*l2*p21^2 - c21*d13*l1*p23^2 - c13*d21*l2*p23^2 - c21*d23*l2*p13^2 - c23*d21*l2*p13^2 + d11*d13*l1^2*p21^2 + d11*d13*l1^2*p23^2 - d21*d23*l2^2*p13^2 - c11*d13*l1*p13*p23 - c13*d11*l1*p13*p23 + c21*d13*l1*p13*p23 - c23*d11*l1*p13*p23 - c11*d23*l2*p13*p23 + c13*d21*l2*p13*p23 + c21*d23*l2*p13*p23 + c23*d21*l2*p13*p23 + d12*l1*p13*p21*p22 - d22*l2*p13*p21*p22 + d11*d23*l1*l2*p13^2 - d13*d21*l1*l2*p21^2 - d13*d21*l1*l2*p23^2 - d11*d13*l1^2*p13*p23 + d21*d23*l2^2*p13*p23 - d11*d23*l1*l2*p13*p23 + d13*d21*l1*l2*p13*p23)/frac;
        s2 = -(c12*c13*p21^2 + c12*c13*p23^2 + c12*c23*p13^2 - c13*c22*p21^2 - c13*c22*p23^2 - c22*c23*p13^2 - c12*c13*p13*p23 - c12*c23*p13*p23 + c13*c22*p13*p23 + c22*c23*p13*p23 - c11*p13*p21*p22 + c21*p13*p21*p22 + c12*d13*l1*p21^2 + c13*d12*l1*p21^2 + c12*d13*l1*p23^2 + c13*d12*l1*p23^2 + c23*d12*l1*p13^2 + c12*d23*l2*p13^2 - c22*d13*l1*p21^2 - c13*d22*l2*p21^2 - c22*d13*l1*p23^2 - c13*d22*l2*p23^2 - c22*d23*l2*p13^2 - c23*d22*l2*p13^2 + d12*d13*l1^2*p21^2 + d12*d13*l1^2*p23^2 - d22*d23*l2^2*p13^2 - c12*d13*l1*p13*p23 - c13*d12*l1*p13*p23 + c22*d13*l1*p13*p23 - c23*d12*l1*p13*p23 - c12*d23*l2*p13*p23 + c13*d22*l2*p13*p23 + c22*d23*l2*p13*p23 + c23*d22*l2*p13*p23 - d11*l1*p13*p21*p22 + d21*l2*p13*p21*p22 + d12*d23*l1*l2*p13^2 - d13*d22*l1*l2*p21^2 - d13*d22*l1*l2*p23^2 - d12*d13*l1^2*p13*p23 + d22*d23*l2^2*p13*p23 - d12*d23*l1*l2*p13*p23 + d13*d22*l1*l2*p13*p23)/frac;
        
        % get T hyp
        t1 = s1*(c13 + d13*l1) - c1*c2*(c11 + d11*l1) - c1*s2*(c12 + d12*l1);
        t2 = s2*(c11 + d11*l1) - c2*(c12 + d12*l1);
        t3 = p13 - c1*(c13 + d13*l1) - c2*s1*(c11 + d11*l1) - s1*s2*(c12 + d12*l1);
        
        if abs(c1*c1 + s1*s1 - 1) < 10e-2 && abs(c2*c2 + s2*s2 - 1) < 10e-2
            sdepth1(i) = l1;
            sdepth2(i) = l2;
            sR(:,(i-1)*3+1:i*3) = [ e11*(c1*c2*f11 - f21*s2 + c2*f31*s1) + e31*(c1*f31 - f11*s1) + e21*(c2*f21 + c1*f11*s2 + f31*s1*s2), e12*(c1*c2*f11 - f21*s2 + c2*f31*s1) + e32*(c1*f31 - f11*s1) + e22*(c2*f21 + c1*f11*s2 + f31*s1*s2), e13*(c1*c2*f11 - f21*s2 + c2*f31*s1) + e33*(c1*f31 - f11*s1) + e23*(c2*f21 + c1*f11*s2 + f31*s1*s2);
                e11*(c1*c2*f12 - f22*s2 + c2*f32*s1) + e31*(c1*f32 - f12*s1) + e21*(c2*f22 + c1*f12*s2 + f32*s1*s2), e12*(c1*c2*f12 - f22*s2 + c2*f32*s1) + e32*(c1*f32 - f12*s1) + e22*(c2*f22 + c1*f12*s2 + f32*s1*s2), e13*(c1*c2*f12 - f22*s2 + c2*f32*s1) + e33*(c1*f32 - f12*s1) + e23*(c2*f22 + c1*f12*s2 + f32*s1*s2);
                e11*(c1*c2*f13 - f23*s2 + c2*f33*s1) + e31*(c1*f33 - f13*s1) + e21*(c2*f23 + c1*f13*s2 + f33*s1*s2), e12*(c1*c2*f13 - f23*s2 + c2*f33*s1) + e32*(c1*f33 - f13*s1) + e22*(c2*f23 + c1*f13*s2 + f33*s1*s2), e13*(c1*c2*f13 - f23*s2 + c2*f33*s1) + e33*(c1*f33 - f13*s1) + e23*(c2*f23 + c1*f13*s2 + f33*s1*s2)];
            sT(:,i) =  [f11*t1 + f21*t2 + f31*t3 - f11*t1_ - f21*t2_ - f31*t3_
                f12*t1 + f22*t2 + f32*t3 - f12*t1_ - f22*t2_ - f32*t3_
                f13*t1 + f23*t2 + f33*t3 - f13*t1_ - f23*t2_ - f33*t3_];
        end
        
    end
    
end
