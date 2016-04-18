% task 1: plot the convergence (amp & pha) vs num of modes 

mappingmodes2D;
clear;close all;clc;
mappingmodes3D;

loadParameters;

Lx = Parameters.Lx;
Ly = Parameters.Ly;
Lz = Parameters.Lz;
Q = Parameters.Q;
Nz = Parameters.Nz;
P = Nz*Q;
% k = Parameters.k; % frequency
rho0 = Parameters.rho0 ;
c0 = Parameters.c0 ;

x0 = 0.1; y0=0.1; z0= -Lz+0.1;
xp = 0.2; yp = 0.3; zp = -Lz +0.4;
xp2 = 1.3;yp2 = 1.4;zp2 = -Lz+1.5;

f = 500;
k = 2*pi*f/340;
Parameters.k = k;
q0 = (10^-4)*rho0*k*c0;
load Modes2D.mat;
load Modes3D.mat;

%------------------------------Forced Response----------------------------%
%--------------------------------H, K, S, M-------------------------------%
K = zeros(P,P); S = zeros(P,1);
for p = 1:P
    u = Modes3D(p,2);
    v = Modes3D(p,3);
    w = Modes3D(p,4);
    ku = u*pi/Lx; kv = v*pi/Ly; kw = w*pi/Lz;
    K(p,p) = ku^2 + kv^2 +kw^2;
    constant = sqrt((2-KronDelta(0,u))/Lx)*sqrt((2-KronDelta(0,v))/Ly)*sqrt((2-KronDelta(0,w))/Lz);
    S(p) = q0*constant*cos(ku*x0)*cos(kv*y0)*cos(kw*z0);
end
H = zeros(P,Q); M = zeros(Q,P);
for p= 1:P
    u = Modes3D(p,2);
    v = Modes3D(p,3);
    w = Modes3D(p,4);
    for q = 1:Q
        m = Modes2D(q,2);
        n = Modes2D(q,3);
        H(p,q) = -1i*k*KronDelta(u,m)*KronDelta(v,n)*sqrt((2-KronDelta(0,w))/Lz);
        M(q,p) = KronDelta(u,m)*KronDelta(v,n)*sqrt((2-KronDelta(0,w))/Lz);
    end
end
Z = getZMat(Parameters);

%--------------------------a----------------------------------------------%
I = eye(P,P);
a = (K + H*(Z\M) -(k^2)*I)\S;
D =K + H*(Z\M);
lambda = sort(eig(D));
V = (Z\M)*a;
% a = (K  -(k^2)*I)\S;
%--------------------------pr_in and pr_out-------------------------------%
pr_in = 0;
for p = 1:P
    u = Modes3D(p,2);
    v = Modes3D(p,3);
    w = Modes3D(p,4);
    ku = u*pi/Lx; kv = v*pi/Ly; kw = w*pi/Lz;
    constant = sqrt((2-KronDelta(0,u))/Lx)*sqrt((2-KronDelta(0,v))/Ly)*sqrt((2-KronDelta(0,w))/Lz);
    pr_in = pr_in + a(p)*constant*cos(ku*xp)*cos(kv*yp)*cos(kw*zp);
end
abs(pr_in)
pr_out = 0;
for q = 1:Q
    m = Modes2D(q,2);
    n = Modes2D(q,3);
    pr_out = pr_out + V(q)*getPrOut_mn(m,n,xp2,yp2,zp2,Parameters);
end
abs(pr_out)
%--------------------bi-orthorgonal modes-(p_in)--------------------------%
phiQ = zeros(Q,1);
for q = 1:Q
    m = Modes2D(q,2);
    n = Modes2D(q,3);    
    phiQ(q) = getPrOut_mn(m,n,xp2,yp2,zp2,Parameters);
end
D =K + H*(Z\M);
[V,lambda_temp] = eig(D);
lambda = zeros(P,1);
for p = 1:P
    lambda(p) = lambda_temp(p,p);
end
exp_mat = [lambda, V.'];
exp_mat = sortrows(exp_mat);
num = 1
for Nn = 1:1:40
% Nn = 25
    pr_recstr_out = 0;
    for n = 1:Nn
        gn = exp_mat(n,2:P+1);gn = gn.';
        lambda_n = exp_mat(n,1);
        c_n = (gn.'*S)/((lambda_n-k^2)*(gn.'*gn));
        pn = phiQ.'*(Z\M)*gn;
        pr_recstr_out = pr_recstr_out + c_n*pn;
    end
    NoM(num) = Nn;
    pr_out_Nn(num) = abs(pr_recstr_out);
    ph_out_Nn(num) = angle(pr_recstr_out);
    num = num +1
end


%--------------------bi-orthorgonal modes-(p_in)--------------------------%
D =K + H*(Z\M);
[V,lambda_temp] = eig(D);
lambda = zeros(P,1);
for p = 1:P
    lambda(p) = lambda_temp(p,p);
end
exp_mat = [lambda, V.'];
exp_mat = sortrows(exp_mat);
zzz = sort(eig(D));
%decide the n^th eigen modes'
% Nn = 15;
% for n = 1:Nn
%     eval(['v',num2str(n),'= exp_mat(',num2str(n),',2:P+1);']); %vn = exp_mat(n,2:P+1);
%     eval(['[maxvalue, index] = max(v',num2str(n),');']); %[maxvalue,index] = max(vn);
%     eval(['u',num2str(n),'= Modes3D(index,2);']); %un = Modes3D(index,2);
%     eval(['v',num2str(n),'= Modes3D(index,2);']); %vn = Modes3D(index,3);
%     eval(['w',num2str(n),'= Modes3D(index,2);']); %wn = Modes3D(index,4);
% end
num = 1
for Nn = 1:1:40
% Nn = 25
    pr_in_recstr = 0;
    for n = 1:Nn
        vn = exp_mat(n,2:P+1);vn = vn.';
        lambda_n = exp_mat(n,1);
        c_n = (vn.'*S)/((lambda_n-k^2)*(vn.'*vn))
        pr_n = 0;
        for p = 1:P
            u = Modes3D(p,2);
            v = Modes3D(p,3);
            w = Modes3D(p,4);
            ku = u*pi/Lx; kv = v*pi/Ly; kw = w*pi/Lz;
            constant = sqrt((2-KronDelta(0,u))/Lx)*sqrt((2-KronDelta(0,v))/Ly)*sqrt((2-KronDelta(0,w))/Lz);
            pr_n = pr_n + vn(p)*constant*cos(ku*xp)*cos(kv*yp)*cos(kw*zp);
        end
        pr_n
        pr_in_recstr = pr_in_recstr + c_n*pr_n;
    end
    NoM(num) = Nn;
    pr_in_Nn(num) = abs(pr_in_recstr);
    ph_in_Nn(num) = angle(pr_in_recstr);
    num = num +1
end


figure(1)
subplot(2,2,1)
plot(NoM,pr_in_Nn) % ,'LineStyle','none','Marker','o'
% xlabel('Number of eigenmodes used for calculation','FontSize',8);
ylabel('Amplitude (pa)','FontSize',8)
grid on;

subplot(2,2,3)
plot(NoM,ph_in_Nn)
xlabel('Number of eigenmodes used for calculation','FontSize',8);
ylabel('Phase','FontSize',8)
grid on;

subplot(2,2,2)
plot(NoM,pr_out_Nn)
ylabel('Amplitude (pa)','FontSize',8)
grid on;
subplot(2,2,4)
plot(NoM,ph_out_Nn)
xlabel('Number of eigenmodes used for calculation','FontSize',8);
ylabel('Phase','FontSize',8)
grid on;
abs(pr_in_recstr)


set(gcf,'PaperUnits','inches','PaperPosition',[0 -1 6 3]);
saveas(gcf, 'cvg_inside.png')
% % % % print('-dtiff','-r300',['processed_figure\',num2str(figure_num),'.tif']);




%------------------------natural property---------------------------------%


