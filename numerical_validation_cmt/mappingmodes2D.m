clear;close all;clc;
loadParameters;
Lx = Parameters.Lx;
Ly = Parameters.Ly;

Nx = 50;
Ny = 50;
Q = Nx*Ny;
Modes2D = zeros(Q,3);
num = 1;
for nx = 0:Nx
    for ny = 0:Ny
            kxyz2 = (nx*pi/Lx)^2 + (ny*pi/Ly)^2;
            Modes2D(num,1) = kxyz2;
            Modes2D(num,2) = nx;
            Modes2D(num,3) = ny;
            num = num+1;
    end
end
Modes2D = sortrows(Modes2D);
save Modes2D.mat Modes2D;

     