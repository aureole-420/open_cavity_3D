function zmat = getZMat( Parameters )

Lx = Parameters.Lx;
Ly = Parameters.Ly;
Lz = Parameters.Lz;
k = Parameters.k; % frequency
rho0 = Parameters.rho0;
c0 = Parameters.c0;
kdir = Parameters.kdir;
Q = Parameters.Q;
% P = Parameters.P;
load Modes2D.mat;
% load Modes3D.mat;
zmat = zeros(Q,Q);
% % % ZMat = zeros(P,Q);
% % % zmat = zeros(Q,Q);
% Paralell computing
for q1 = 1:Q
    for q = 1:Q
        m = Modes2D(q,2);
        n = Modes2D(q,3);
        m1 = Modes2D(q1,2);
        n1 = Modes2D(q1,3);
%         if ~exist('zmat','dir')
%             mkdir('zmat');
        if ~exist(kdir,'dir')
            mkdir(kdir);
        end
%         a = exist(['zmat/z_',num2str(m),'_',num2str(n),'_',num2str(m1),'_',num2str(n1),'.mat'],'file');
%         b = exist(['zmat/z_',num2str(m1),'_',num2str(n1),'_',num2str(m),'_',num2str(n),'.mat'],'file');
        a = exist([kdir,'/z_',num2str(m),'_',num2str(n),'_',num2str(m1),'_',num2str(n1),'.mat'],'file');
        b = exist([kdir,'/z_',num2str(m1),'_',num2str(n1),'_',num2str(m),'_',num2str(n),'.mat'],'file');

        if a
            load([kdir,'/z_',num2str(m),'_',num2str(n),'_',num2str(m1),'_',num2str(n1),'.mat']);
            eval(['zmat(q1,q)=z_',num2str(m),'_',num2str(n),'_',num2str(m1),'_',num2str(n1),';']);
%             load(['zmat/z_',num2str(m),'_',num2str(n),'_',num2str(m1),'_',num2str(n1),'.mat']);
%             eval(['zmat(q1,q)=z_',num2str(m),'_',num2str(n),'_',num2str(m1),'_',num2str(n1),';']);
       
        elseif b
            load([kdir,'/z_',num2str(m1),'_',num2str(n1),'_',num2str(m),'_',num2str(n),'.mat']);
            eval(['zmat(q1,q)=z_',num2str(m1),'_',num2str(n1),'_',num2str(m),'_',num2str(n),';']); 
%             load(['zmat/z_',num2str(m1),'_',num2str(n1),'_',num2str(m),'_',num2str(n),'.mat']);
%             eval(['zmat(q1,q)=z_',num2str(m1),'_',num2str(n1),'_',num2str(m),'_',num2str(n),';']); 
        
        else
            zmat(q1,q) = getZmn(m,m1,n,n1,Parameters);
            eval(['z_',num2str(m),'_',num2str(n),'_',num2str(m1),'_',num2str(n1),'=zmat(q1,q);'])
%             save zmat/z_m_n_m1_n1.mat z_m_n_m1_n1;
            save([kdir,'/z_',num2str(m),'_',num2str(n),'_',num2str(m1),'_',num2str(n1),'.mat '],['z_',num2str(m),'_',num2str(n),'_',num2str(m1),'_',num2str(n1)]);
%             save(['zmat/z_',num2str(m),'_',num2str(n),'_',num2str(m1),'_',num2str(n1),'.mat '],['z_',num2str(m),'_',num2str(n),'_',num2str(m1),'_',num2str(n1)]);
        end
    end
end

end

