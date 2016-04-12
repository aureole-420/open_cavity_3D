loadParameters;
% f = 500;
% k = 2*pi*f/340;
Lx = Parameters.Lx;
Ly = Parameters.Ly;
a = sqrt(Lx*Ly);
k = 0.1/a;
Parameters.k = k;
m = 0;
m1 = 0;
n = 0;
n1 =0;
getZmn(m,m1,n,n1,Parameters)