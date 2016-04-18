function [ result ] = getPrOut_mn(m,n,xp,yp,zp,Parameters)
%GETPROUT Summary of this function goes here
%   Detailed explanation goes here
Lx = Parameters.Lx;
Ly = Parameters.Ly;
k = Parameters.k;
rho0 = Parameters.rho0;
c0 = Parameters.c0;
Q = Parameters.Q;

func = @(x1,y1) pr_out_mn_int(x1,y1,m,n,xp,yp,zp,Parameters);

result = integral2(func,0,Lx,0,Ly,'AbsTol',0,'RelTol',10^-8);


end

