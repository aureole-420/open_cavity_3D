function [ Zmn ] = getZmn(m,m1,n,n1,Parameters)
% Return the radiation impedance Zmn relating
% (m,n) mode and (m1,n1) mode
Lx = Parameters.Lx;
Ly = Parameters.Ly;
k = Parameters.k;
rho0 = Parameters.rho0;
c0 = Parameters.c0;

func = @(x,y) gxmm1(m,m1,x,Parameters).*gynn1(n,n1,y,Parameters).*exp(1i*k*sqrt(x.^2+y.^2))./(sqrt(x.^2+y.^2));
% func = @(x,y) (-1i*k/(2*pi*Lx*Ly))*func0(x,y);
constant = (-1i*k/(2*pi*Lx*Ly))*sqrt(2-KronDelta(0,m))*sqrt(2-KronDelta(0,m1));
constant = constant*sqrt(2-KronDelta(0,n))*sqrt(2-KronDelta(0,n1));
Zmn = constant*integral2(func,0,Lx,0,Ly,'AbsTol',0,'RelTol',10^-8);

end

