function [] = plotXY(gn, z_plot, Parameters)
% plot pressure distribution of XY plane at z= z0
Lx = Parameters.Lx;
Ly = Parameters.Ly;
Lz = Parameters.Lz;
k = Parameters.k; % frequency

Q = Parameters.Q;
Nz = Parameters.Nz; % Nz
P = Q*Nz;

numx = 500;
numy = 500;
x = 0:Lx/numx:Lx;
y = 0:Ly/numy:Ly;
pressure_field = zeros(length(y),length(x));

load Modes3D.mat;
for p = 1:P
    u = Modes3D(p,2);
    v = Modes3D(p,3);
    w = Modes3D(p,4);
    ku = u*pi/Lx;kv = v*pi/Ly;kw = w*pi/Lz;
    constant = sqrt((2-KronDelta(0,u))/Lx)*sqrt((2-KronDelta(0,v))/Ly)*sqrt((2-KronDelta(0,w))/Lz);
    pressure_field  = pressure_field + gn(p)*(constant*cos(kw*z_plot)*cos(ku*(x.'))*cos(kv*y)).';
end
mesh(x,y,abs(pressure_field));
xlim([0,Lx]);
ylim([0,Ly]);
axis equal
xlabel('x (m)');ylabel('y (m)');
colorbar;
colormap jet;
view(2);
saveas(gcf,'XY.png')

end

