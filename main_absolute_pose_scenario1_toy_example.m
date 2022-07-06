clc;clear all;close all;

addpath utils/
addpath solvers\
addpath libso3\

K = [4800 0 960;0 4800 540;0 0 1];

simu_noise_level = 0.1;

num_features = 100;% always use 100 features for pose estimation

% 7 methods for regular case
methods = {'gP+s','MMPPE','AGW','gPnP','5pt','PL','PL+O'};%
 
% random draw a camera optical axis
opt_z = [rand(1)-0.5,rand(1)-0.5,-1]';
opt_z = opt_z ./ norm(opt_z);
% find real rotation
rotaxis = cross([0;0;1],opt_z);
rotaxis = rotaxis ./ norm(rotaxis);
theta = acos(dot(opt_z,[0;0;1]));
R = expSO3(theta*rotaxis);
% verify the real rotation is correct
assert(abs(dot(R*[0;0;1],opt_z)-1)<1e-6);

% camera pose: ground truth, since the camera should be above the plane
R1 = R; t1 = [0 0 1+rand(1)*0.5]';
gtR = R1;gtt = t1;

n1 = [0;0;1];d1 = 0;
n1 = n1./norm(n1);

refra_n1 = 1.0;% air refractive index;
refra_n2 = 1.5;% simulated refractive index;
miu = refra_n2/refra_n1;%

[q,Q,ref_As,planes,tr1s]=dataGenerator(num_features,R1,n1,t1,refra_n1,refra_n2,simu_noise_level,K);
            
figure;colormap jet;
plot3(t1(1),t1(2),t1(3),'k-d','MarkerSize',15);hold on;
h1=plotCameraFrustum(R1,-R1*t1,K,0.5,'r',1.5);
for ii = 1:num_features
    h=line([t1(1) ; ref_As(1,ii)], [t1(2) ; ref_As(2,ii)],[t1(3) ; ref_As(3,ii)] ,'LineWidth', 2 , 'Color', 'b');hold on;
    h.Color(4) = 0.5;
    h2=plot3(ref_As(1,ii),ref_As(2,ii),ref_As(3,ii),'y-s','MarkerFaceColor','y');hold on;
    h=line([Q(1,ii) ; ref_As(1,ii)], [Q(2,ii) ; ref_As(2,ii)],[Q(3,ii) ; ref_As(3,ii)] ,'LineWidth', 2, 'Color', 'g');
    h.Color(4) = 0.5;
end
[xx,yy] = meshgrid(min(ref_As(1,:))-1:0.1:max(ref_As(1,:))+1,min(ref_As(2,:))-1:0.1:max(ref_As(1,:))+1);
zz = (-d1 - n1(1)*xx - n1(2)*yy)/n1(3);
h3=mesh(xx,yy,zz);hold on
h4=plot3(Q(1,:),Q(2,:),Q(3,:),'ko','MarkerFaceColor','k');hold on;
xlabel('x: (m)','Interpreter','latex');ylabel('y: (m)','Interpreter','latex');zlabel('z: (m)','Interpreter','latex');
title('Simulation case: regular','Interpreter','latex');
legend([h1,h2,h3,h4],{'Camera','Refractive Point','Refractive Plane','World Point'},'Interpreter','latex');

%
n1c = R1 * -n1;
d = abs(dot(n1,t1));

for kk = 1:length(methods)
    if strcmp(methods{kk},'gP+s')
        addpath './thirdparty/genposeandscale-master'
        try
            tic;
            [Rana,tana] = gPplus_Refraction(Q,q,K,n1c,d,refra_n1,refra_n2,R1,t1);
            duration = toc;
        catch
            Rana = eye(3);
            tana = zeros(3,1);
        end
        disp(['ratation error of gP+s: ', num2str(norm(Rana'*R1-eye(3)))]);
        disp(['translation error of gP+s: ', num2str(norm(-Rana'*tana-t1))]);
        rmpath('./thirdparty/genposeandscale-master');
    end

    if strcmp(methods{kk},'MMPPE')
        addpath(genpath('./thirdparty/MinimalMultiPerspectivePose-master'));
        try
            tic;
            [Rana,tana] = gpnp_eccv18_Refraction(Q,q,K,n1c,d,refra_n1,refra_n2,R1,t1);
            duration = toc;
        catch
            Rana = eye(3);
            tana = zeros(3,1);
        end
        disp(['ratation error of MMPPE: ', num2str(norm(Rana'*R1-eye(3)))]);
        disp(['translation error of MMPPE: ', num2str(norm(-Rana'*tana-t1))]);
        rmpath(genpath('./thirdparty/MinimalMultiPerspectivePose-master'));
    end 


    if strcmp(methods{kk},'AGW')
        addpath(genpath('./thirdparty/AgrawalCVPR12WebCode/RealDataWebCode'));
        try
            tic;
            [Ragw,Cagw] = abs_pose_agw(Q,q,K,refra_n1,refra_n2,1,0,0);
            duration = toc;
        catch
            Ragw = eye(3);
            Cagw = zeros(3,1);
        end
        disp(['ratation error of AGW: ', num2str(norm(Ragw'*R1-eye(3)))]);
        disp(['translation error of AGW: ', num2str(norm(Cagw-t1))]);
        rmpath(genpath('./thirdparty/AgrawalCVPR12WebCode/RealDataWebCode'));
    end 

    if strcmp(methods{kk},'gPnP')
        try
            tic
            [Rana,tana] = gPnP_Refraction(Q,q,K,n1c,d,refra_n1,refra_n2);
            duration = toc;
        catch
            disp(['Please install opengv!']);
            Rana = eye(3);
            tana = zeros(3,1);
        end
        disp(['ratation error of gPnP: ', num2str(norm(Rana'*R1-eye(3)))]);
        disp(['translation error of gPnP: ', num2str(norm(-Rana'*tana-t1))]);
    end

    if strcmp(methods{kk},'uPnP')
        try
            tic
            [Rana,tana] = uPnP_Refraction(Q,q,K,n1c,d,refra_n1,refra_n2);
            duration = toc;
        catch
            disp(['Please install opengv!'])
            Rana = eye(3);
            tana = zeros(3,1);
        end
        disp(['ratation error of uPnP: ', num2str(norm(Rana'*R1-eye(3)))]);
        disp(['translation error of uPnP: ', num2str(norm(-Rana'*tana-t1))]);
    end

    if strcmp(methods{kk},'5pt')
        addpath('thirdparty/refractive_pose-master')
        try
            tic
            [Rana,Cana] = Hanner_8pt(Q,q,K,miu,R1,t1);
            duration = toc;
        catch
            Rana = eye(3);
            Cana = zeros(3,1);
        end
        rmpath('thirdparty/refractive_pose-master')
        disp(['ratation error of 5pt: ', num2str(norm(Rana'*R1-eye(3)))]);
        disp(['translation error of 5pt: ', num2str(norm(Cana-t1))]);
    end

    if strcmp(methods{kk},'PL')
        try
            tic
            [Rana,tana] = PnP_Refraction(Q,q,K,n1c,d,refra_n1,refra_n2);
            duration = toc;
        catch
            Rana = eye(3);
            tana = zeros(3,1);
        end
        disp(['ratation error of PL: ', num2str(norm(Rana'*R1-eye(3)))]);
        disp(['translation error of PL: ', num2str(norm(-Rana'*tana-t1))]);
    end

    if strcmp(methods{kk},'PL+O')
        try
            tic
            [Rana,tana] = PnP_Refraction(Q,q,K,n1c,d,refra_n1,refra_n2,1);
            duration = toc;
        catch
            Rana = eye(3);
            tana = zeros(3,1);
        end
        disp(['ratation error of PL+O: ', num2str(norm(Rana'*R1-eye(3)))]);
        disp(['translation error of PL+O: ', num2str(norm(-Rana'*tana-t1))]);
    end
end
