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
NUM = 1
for k = 9.24:0.05:9.24
    ksource(NUM) = k;
    freq(NUM) = k*340/2/pi;
    Parameters.k = k;
    Parameters.kdir = ['zmat/k',num2str(k)];
    q0 = (10^-4)*rho0*k*c0;
%     q0 = 4*pi*10^-4;
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

    %--------------------------pr_in and pr_out-------------------------------%
% % %     pr_in = 0;
% % %     for p = 1:P
% % %         u = Modes3D(p,2);
% % %         v = Modes3D(p,3);
% % %         w = Modes3D(p,4);
% % %         ku = u*pi/Lx; kv = v*pi/Ly; kw = w*pi/Lz;
% % %         constant = sqrt((2-KronDelta(0,u))/Lx)*sqrt((2-KronDelta(0,v))/Ly)*sqrt((2-KronDelta(0,w))/Lz);
% % %         pr_in = pr_in + a(p)*constant*cos(ku*xp)*cos(kv*yp)*cos(kw*zp);
% % %     end
% % %     abs(pr_in)
% % %     pr_out = 0;
% % %     for q = 1:Q
% % %         m = Modes2D(q,2);
% % %         n = Modes2D(q,3);
% % %         pr_out = pr_out + V(q)*getPrOut_mn(m,n,xp2,yp2,zp2,Parameters);
% % %     end
% % %     abs(pr_out)
% % %     
% % %     PR_OUT(NUM) = pr_out;
% % %     PR_IN(NUM) = pr_in;
       
% % % % %     % -----------------bi-ortho modes (using Nm =20 eigen-modes)--------------%
% % % % %     
% % % % %     
% % % % %     % % % % D =K + H*(Z\M);
% % % % %     [V,lambda_temp] = eig(D);
% % % % %     lambda = zeros(P,1);
% % % % %     for p = 1:P
% % % % %         lambda(p) = lambda_temp(p,p);
% % % % %     end
% % % % %     exp_mat = [lambda, V.'];
% % % % %     exp_mat = sortrows(exp_mat);
% % % % %     
% % % % %     phiQ = zeros(Q,1);
% % % % %     for q = 1:Q
% % % % %         m = Modes2D(q,2);
% % % % %         n = Modes2D(q,3);
% % % % %         phiQ(q) = getPrOut_mn(m,n,xp2,yp2,zp2,Parameters);
% % % % %     end
% % % % %     
% % % % %     Nn = 20;
% % % % %     pr_in_recstr = 0;
% % % % %     pr_out_recstr = 0;
% % % % %     for n = 1:Nn
% % % % %         gn = exp_mat(n,2:P+1);gn = gn.';
% % % % %         lambda_n = exp_mat(n,1);
% % % % %         c_n = (gn.'*S)/((lambda_n-k^2)*(gn.'*gn));
% % % % %         % reconstructed pr_out
% % % % %         pn = phiQ.'*(Z\M)*gn;
% % % % %         pr_out_recstr = pr_out_recstr + c_n*pn;
% % % % %         pr_n = 0;
% % % % %         pn = 0;
% % % % %         % reconstructed pr_in
% % % % %         for p = 1:P
% % % % %             u = Modes3D(p,2);
% % % % %             v = Modes3D(p,3);
% % % % %             w = Modes3D(p,4);
% % % % %             ku = u*pi/Lx; kv = v*pi/Ly; kw = w*pi/Lz;
% % % % %             constant = sqrt((2-KronDelta(0,u))/Lx)*sqrt((2-KronDelta(0,v))/Ly)*sqrt((2-KronDelta(0,w))/Lz);
% % % % %             pr_n = pr_n + gn(p)*constant*cos(ku*xp)*cos(kv*yp)*cos(kw*zp);
% % % % %         end
% % % % %         pr_in_recstr = pr_in_recstr + c_n*pr_n;
% % % % %     end
% % % % %     pr_in_spl(NUM) = 20*log10(abs(pr_in_recstr)/(2*10^-5));
% % % % %     pr_out_spl(NUM) = 20*log10(abs(pr_out_recstr)/(2*10^-5));
    
    
    
    
    
    NUM=NUM+1
end
% % % % 
% % % % load freqfea.mat;load pin_spl.mat;load pout_spl.mat;
% % % % figure(1)
% % % % % plot(freqfea,pin_spl,'r','LineWidth',1,'LineStyle','-')
% % % % plot(freqfea,pin_spl,'r','LineWidth',2)
% % % % hold on
% % % % % plot(freq, pr_in_spl,'b','LineWidth',0.5,'LineStyle','none','Marker','o','MarkerSize',5)
% % % % plot(freq, pr_in_spl,'b','LineWidth',1)
% % % % hold on;
% % % % plot(freqfea,pout_spl,'k','LineWidth',2,'LineStyle','-')
% % % % % plot(freqfea,pout_spl,'g','LineWidth',1,'LineStyle','-')
% % % % hold on;
% % % % plot(freq,pr_out_spl,'g','LineWidth',1)
% % % % % plot(freq,pr_out_spl,'b','LineWidth',0.5,'LineStyle','none','Marker','x','MarkerSize',5)
% % % % xlim([30 500]);
% % % % ylim([-20 60])
% % % % xlabel('Frequency (Hz)','FontSize',8);
% % % % ylabel('SPL (dB)','FontSize',8);
% % % % h_legend = legend('COMSOL (inside)','MATLAB (inside)','COMSOL (outside)','MATLAB (outside)','Location','SouthWest');
% % % % set(h_legend,'FontSize',8)
% % % % set(gcf,'PaperUnits','inches','PaperPosition',[0 -1 4 3]);
% % % % saveas(gcf, 'prediction_20eigenmodes.png')

% % % % using P eigen modes
% % % PR_OUT_spl = 20*log10(abs(PR_OUT)/(2*10^-5))
% % % PR_IN_spl = 20*log10(abs(PR_IN)/(2*10^-5))
% % % plot(freq,PR_IN_spl)
% % % hold on;
% % % plot(freq,PR_OUT_spl)
% % % legend('in','out')




% % % % %-----------------bi-orthorgonal modes-(p_in & p_out)---------------------%
D =K + H*(Z\M);
aaa = sort(sqrt(eig(D)))*340/2/pi;
[V,lambda_temp] = eig(D);
lambda = zeros(P,1);
for p = 1:P
    lambda(p) = lambda_temp(p,p);
end
exp_mat = [lambda, V.'];
exp_mat = sortrows(exp_mat);

phiQ = zeros(Q,1);
for q = 1:Q
    m = Modes2D(q,2);
    n = Modes2D(q,3);
    phiQ(q) = getPrOut_mn(m,n,xp2,yp2,zp2,Parameters);
end
% % % % 
% ---------- plot  modal shape for nnn = 6th eigenmode------------%
% % % gn = exp_mat(1,2:P+1);gn = gn.';
% % % plotXY(gn,0,Parameters)
% % % plotXZ(gn,0,Parameters)
% % % plotYZ(gn,0,Parameters)
% plotXY(a,-Lz,Parameters)
% % plotXZ(a,Ly,Parameters)
% % plotYZ(a,0,Parameters)


% % % % num = 1
% % % % %%%---------------------------(1)plot Pr~Nn-------------------------------%
% % % % %%%%-------------------------(2)plot cn~n---------------------------------%
% % % % % bi_basis
% % % % % orth_basis
% % % % for Nn = 1:1:40
% % % % % Nn = 25
% % % %     pr_in_recstr = 0;
% % % %     pr_out_recstr = 0;
% % % %     for n = 1:Nn
% % % %         gn = exp_mat(n,2:P+1);gn = gn.';
% % % %         lambda_n = exp_mat(n,1);
% % % %         c_n = (gn.'*S)/((lambda_n-k^2)*(gn.'*gn));
% % % %         bi_basis(n) = c_n;
% % % %         orth_basis(n) = a(n);
% % % %         % reconstructed pr_out
% % % %         pn = phiQ.'*(Z\M)*gn;
% % % %         pr_out_recstr = pr_out_recstr + c_n*pn;
% % % %         pr_n = 0;
% % % %         pn = 0;
% % % %         % reconstructed pr_in
% % % %         for p = 1:P
% % % %             u = Modes3D(p,2);
% % % %             v = Modes3D(p,3);
% % % %             w = Modes3D(p,4);
% % % %             ku = u*pi/Lx; kv = v*pi/Ly; kw = w*pi/Lz;
% % % %             constant = sqrt((2-KronDelta(0,u))/Lx)*sqrt((2-KronDelta(0,v))/Ly)*sqrt((2-KronDelta(0,w))/Lz);
% % % %             pr_n = pr_n + gn(p)*constant*cos(ku*xp)*cos(kv*yp)*cos(kw*zp);
% % % %         end
% % % %         pr_n
% % % %         pr_in_recstr = pr_in_recstr + c_n*pr_n;
% % % %     end
% % % %     NoM(num) = Nn;
% % % %     pr_in_Nn(num) = abs(pr_in_recstr);
% % % %     ph_in_Nn(num) = angle(pr_in_recstr);
% % % %     pr_out_Nn(num) = abs(pr_out_recstr);
% % % %     ph_out_Nn(num) = angle(pr_out_recstr);
% % % %     num = num +1
% % % % end
% % % % figure(1)
% % % % nmodes = 1:40;
% % % % % plot(nmodes, abs(bi_basis),'LineWidth',2)
% % % % % hold on;
% % % % % plot(nmodes, abs(orth_basis),'LineWidth',2);
% % % % bar1 = bar(nmodes,[abs(orth_basis);abs(bi_basis)]','BarWidth',1);
% % % % set(bar1(1),'FaceColor',[0 0 0]);
% % % % set(bar1(2),'FaceColor',[1 0 0]);
% % % % legend('Orthogonal bases','Bi-orthogonal bases')
% % % % xlabel('Order of Modes');
% % % % ylabel('Modal coefficients')
% % % % xlim([0,30])
% % % % set(gcf,'PaperUnits','inches','PaperPosition',[0 -1 4 3]);
% % % % saveas(gcf, 'modal_coefficients.png')

% % % % % % figure(1)
% % % % % % subplot(2,2,1)
% % % % % % plot(NoM,pr_in_Nn,'LineWidth',2) % ,'LineStyle','none','Marker','o'
% % % % % % % xlabel('Number of eigenmodes used for calculation','FontSize',8);
% % % % % % ylabel('Amplitude (pa)','FontSize',8)
% % % % % % grid on;
% % % % % % 
% % % % % % subplot(2,2,3)
% % % % % % plot(NoM,ph_in_Nn,'LineWidth',2)
% % % % % % xlabel('Number of eigenmodes used for calculation','FontSize',8);
% % % % % % ylabel('Phase','FontSize',8)
% % % % % % grid on;
% % % % % % 
% % % % % % subplot(2,2,2)
% % % % % % plot(NoM,pr_out_Nn,'LineWidth',2)
% % % % % % ylabel('Amplitude (pa)','FontSize',8)
% % % % % % grid on;
% % % % % % subplot(2,2,4)
% % % % % % plot(NoM,ph_out_Nn,'LineWidth',2)
% % % % % % xlabel('Number of eigenmodes used for calculation','FontSize',8);
% % % % % % ylabel('Phase','FontSize',8)
% % % % % % grid on;
% % % % % % abs(pr_in_recstr)
% % % % % % 
% % % % % % 
% % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0 -1 6 3]);
% % % % % % saveas(gcf, 'cvg_inside.png')
% % % % % % % print('-dtiff','-r300',['processed_figure\',num2str(figure_num),'.tif']);


%------------------------natural property---------------------------------%


