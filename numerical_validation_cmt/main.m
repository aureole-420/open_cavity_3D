
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
xp = 0.2; yp = 0.3; zp = -Lz +0.3;

f = 500;
k = 2*pi*f/340;
Parameters.k = k;
q0 = (10^-4)*rho0*k*c0;
load Modes2D.mat;
load Modes3D.mat;

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
% a = (K  -(k^2)*I)\S;    
%--------------------------pr---------------------------------------------%
pr = 0;
for p = 1:P
    u = Modes3D(p,2);
    v = Modes3D(p,3);
    w = Modes3D(p,4);
    ku = u*pi/Lx; kv = v*pi/Ly; kw = w*pi/Lz;
    constant = sqrt((2-KronDelta(0,u))/Lx)*sqrt((2-KronDelta(0,v))/Ly)*sqrt((2-KronDelta(0,w))/Lz);
    pr = pr + a(p)*constant*cos(ku*xp)*cos(kv*yp)*cos(kw*zp);
end
abs(pr)



% % % 
% % % %------------------------Zmat, S, P --------------------------------------%
% % % Zmat = getZMat(Parameters);
% % % S = zeros(Q,1);P = zeros(Q,Q);
% % % for q = 1:Q
% % %     u = Modes2D(q,2);
% % %     v = Modes2D(q,3);
% % %     ku = u*pi/Lx; kv = v*pi/Ly;
% % %     tempS = 0;tempP = 0;
% % %     for w = 0:Nz-1
% % %          kw = w*pi/Lz;
% % %         constant =sqrt((2-KronDelta(0,u))/Lx)*sqrt((2-KronDelta(0,v))/Ly)*((2-KronDelta(0,w))/Lz);
% % %         tempS = tempS + q0*constant*cos(ku*x0)*cos(kv*y0)*cos(kw*z0)/(ku^2+kv^2+kw^2 -k^2);
% % %         tempP = tempP -1i*k*((2-KronDelta(0,w))/Lz)/(ku^2+kv^2+kw^2 -k^2);
% % %     end
% % %     P(q,q) = tempP;
% % %     S(q) = tempS;
% % % end
% % % V = (P+Zmat)\S;
% % % 
% % % %-------------------------pa(xp,yp,zp)------------------------------------%
% % % pr = 0;
% % % for q = 1:Q
% % %     u = Modes2D(q,2);
% % %     v = Modes2D(q,3);
% % %     temp1 = 0; temp2 = 0;
% % %     for w = 0:Nz-1
% % %         ku = u*pi/Lx; kv = v*pi/Ly; kw = w*pi/Lz;
% % %         constant1 = ((2-KronDelta(0,u))/Lx)*((2-KronDelta(0,v))/Ly)*((2-KronDelta(0,w))/Lz);
% % %         constant2 = sqrt((2-KronDelta(0,u))/Lx)*sqrt((2-KronDelta(0,v))/Ly)*((2-KronDelta(0,w))/Lz);
% % %         temp1 = temp1 + q0*constant1*cos(ku*xp)*cos(kv*yp)*cos(kw*zp)*cos(ku*x0)*cos(kv*y0)*cos(kw*z0)/(ku^2+kv^2+kw^2-k^2);
% % %         temp2 = temp2 +1i*k*V(q)*constant2*cos(ku*xp)*cos(kv*yp)*cos(kw*zp)/(ku^2+kv^2+kw^2 - k^2);
% % %     end
% % %     pr = pr +temp1 +temp2;
% % % end
% % % abs(pr)
% % % VV = abs(V);

