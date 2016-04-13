clear;close all;clc;
loadParameters;
Lx = Parameters.Lx;
Ly = Parameters.Ly;
Lz = Parameters.Lz;

load Modes2D.mat;
Q = Parameters.Q;
Nz = Parameters.Nz;
Modes3D = zeros(Q*Nz,4);
num = 1;
for q = 1:Q
    nx = Modes2D(q,2);
    ny = Modes2D(q,3);
        for nz = 0:Nz-1
            kxyz2 = (nx*pi/Lx)^2 + (ny*pi/Ly)^2 +(nz*pi/Lz)^2;
            Modes3D(num,1) = kxyz2;
            Modes3D(num,2) = nx;
            Modes3D(num,3) = ny;
            Modes3D(num,4) = nz;
            num = num+1;
    end
end
Modes3D = sortrows(Modes3D);
save Modes3D.mat Modes3D;

     