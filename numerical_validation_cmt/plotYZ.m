function [] = plotYZ(gn, x_plot, Parameters)
% plot pressure distribution of XY plane at z= z0
Lx = Parameters.Lx;
Ly = Parameters.Ly;
Lz = Parameters.Lz;

k = Parameters.k; % frequency

Q = Parameters.Q;
Nz = Parameters.Nz; % Nz
P = Q*Nz;


numz = 500;
numy = 500;
z = -Lz:Lz/numz:-0;
y = 0:Ly/numy:Ly;
pressure_field = zeros(length(z),length(y));
load Modes3D.mat;
for p = 1:P
    u = Modes3D(p,2);
    v = Modes3D(p,3);
    w = Modes3D(p,4);
    ku = u*pi/Lx;kv = v*pi/Ly;kw = w*pi/Lz;
    constant = sqrt((2-KronDelta(0,u))/Lx)*sqrt((2-KronDelta(0,v))/Ly)*sqrt((2-KronDelta(0,w))/Lz);
    pressure_field  = pressure_field + gn(p)*(constant*cos(ku*x_plot)*cos(kv*(y.'))*cos(kw*z)).';
end
mesh(y,z,abs(pressure_field));
xlim([0,Ly]);
ylim([-Lz,0]);
xlabel('y (m)');ylabel('z (m)');
axis equal;
colorbar;
colormap jet;
view(2);
saveas(gcf,'YZ.png')
end

