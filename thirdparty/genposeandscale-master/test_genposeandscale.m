N = 20;
%%
% generate random camera centers and points
p = rand(3,N);
q = 10*rand(3,N);

% get rays between cameras and points
d = q-p;
d = d./repmat(sqrt(sum(d.^2)),3,1);

%%
% generate random rotation
R = rand(3);
[U,~,V] = svd(R);
R = U*diag([1,1,1])*V';
if det(R)<0,
    R = -R;
end

% generate random translation
t = rand(3,1);

% generate random scale
s = rand;

%%
% apply pose and scale
Pq = R*q+repmat(t,1,N);
sp = s*p;

% [R1 t1] = gPPnP(d',sp',Pq',1e-6,1);
% R1 = R1';
% t1 = -R1*t1;

%%
% use transformed input in solver
[Psolns,ssolns] = genposeandscale(sp,d,Pq);

%%
% check answer
fprintf('---- results ----\n');

fprintf('Ground truth pose:\n');    
disp([R t]);

fprintf('Ground truth scale:\n');    
disp(s);

nsolns = length(Psolns);
for i=1:nsolns,
    Rsoln = Psolns{i}(1:3,1:3);
    tsoln = Psolns{i}(1:3,4);
    ssoln = ssolns(i);
    
    fprintf('Pose solution %d:\n',i);
    disp([Rsoln tsoln]); 
    
    fprintf('Scale solution %d:\n',i);
    disp(ssoln);
end