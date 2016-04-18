function [ result ] = pr_out_mn_int(x1,y1,m,n,xp,yp,zp,Parameters)

Lx = Parameters.Lx;
Ly = Parameters.Ly;
k = Parameters.k;
rho0 = Parameters.rho0;
c0 = Parameters.c0;

constant = -1i*k*sqrt((2-KronDelta(0,m))/Lx)*sqrt((2-KronDelta(0,n))/Ly);

r = sqrt((xp-x1).^2+(yp-y1).^2+(zp^2));
km = m*pi/Lx;
kn = n*pi/Ly;
result = constant*cos(km*x1).*cos(kn*y1).*exp(1i*k*r)./(2*pi*r);

end

