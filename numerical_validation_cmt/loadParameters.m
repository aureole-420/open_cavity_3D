% Geometrical parameters: open cavity
% Parameters.Lx = 0.3;
% Parameters.Ly = 0.5;
% Parameters.Lz = 0.8;

Parameters.Lx = 0.432;
Parameters.Ly = 0.67;
Parameters.Lz = 0.598;

% Parameters.Lx = 1;
% Parameters.Ly = 1;
% Parameters.Lz = 1;

% % Parameters.Nx = 20; % discretized
% % Parameters.Ny = 30; % Ny
Parameters.Nz = 10; % Nz

% Parameters.P = 200; % modes for the cavity
Parameters.Q = 40; % radiation modes for the open space

Parameters.k = 10*2*pi/340; % frequency

Parameters.rho0 = 1.225;
Parameters.c0 = 340;
Parameters.kdir = ['zmat/k',num2str(Parameters.k)];
