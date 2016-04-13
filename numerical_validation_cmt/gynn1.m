function [ r ] = gynn1(n,n1,y,Parameters )
Ly = Parameters.Ly;
kn = n*pi/Ly;
kn1 = n1*pi/Ly;

if n == n1 && n == 0
    r = 2*(Ly - y);
elseif n==n1 && n~=0
    r = - (1/kn)*sin(kn.*y) + (Ly-y).*cos(kn*y);
elseif n ~= n1
    r = epsilon(n,n1)*(kn1*sin(kn1*y) - kn*(sin(kn*y)))./(kn^2 -kn1^2);
end


end

