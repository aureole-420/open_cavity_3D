function [ r ] = gxmm1(m,m1,x,Parameters )
Lx = Parameters.Lx;
km = m*pi/Lx;
km1 = m1*pi/Lx;

if m == m1 && m == 0
    r = 2*(Lx - x);
elseif m==m1 && m~=0
    r = - (1/km)*sin(km.*x) + (Lx-x).*cos(km*x);
elseif m ~= m1
    r = epsilon(m,m1)*(km1*sin(km1*x) - km*(sin(km*x)))./(km^2 -km1^2);
end


end

