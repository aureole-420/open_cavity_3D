function [] = plotXZ(gn, y_plot, Parameters)
% plot pressure distribution of XY plane at z= z0
Lx = Parameters.Lx;
Ly = Parameters.Ly;
Lz = Parameters.Lz;
k = Parameters.k; % frequency

Q = Parameters.Q;
Nz = Parameters.Nz; % Nz
P = Q*Nz;

numx = 500;
numz = 500;
x = 0:Lx/numx:Lx;
z = -Lz:Lz/numz:0;
pressure_field = zeros(length(z),length(x));
load Modes3D.mat;
for p = 1:P
    u = Modes3D(p,2);
    v = Modes3D(p,3);
    w = Modes3D(p,4);
    ku = u*pi/Lx;kv = v*pi/Ly;kw = w*pi/Lz;
    constant = sqrt((2-KronDelta(0,u))/Lx)*sqrt((2-KronDelta(0,v))/Ly)*sqrt((2-KronDelta(0,w))/Lz);
    pressure_field  = pressure_field + gn(p)*(constant*cos(kv*y_plot)*cos(ku*(x.'))*cos(kw*z)).';
end

mesh(x,z,abs(pressure_field));
xlim([0,Lx]);
ylim([-Lz,0]);
xlabel('x (m)');ylabel('z (m)');
axis equal;
colorbar;
colormap jet;
view(2);
saveas(gcf,'XZ.png')
end

