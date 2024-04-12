%% Parameters
N =30;dx=0.05;%dx=0.0075707;
x = (-N/2+1:N/2);
x = dx*x;

%Physical Constants
M = 1;                              %Mass used for calculation
m_e = 1;                            %Electron Mass
hbar = 1;                           %Planck
N_f = 2;                            %Flavour
N_c = 3;                            %Ferm Const
gamma_E = 1;                        %
alpha_s = 1;                        %Strong coupling
g = 1;                              %
Temp = 1;                           %Temperature
C_f = 1;                            %
m_D = g^2*Temp^2/3 *(N_c+N_f/2);    %Quark Mass
a_0 = 1;                            %Bohr Radius

%% Matrices and non sparse (1D) Calculation

%Potentials
% V = -diag(1./x); % Coulomb
% %V = 0.5*diag(x.^2); %Oscillator
% %V = diag(x.*0); % Box
% 
% % %%Quark Antiquark
% %M = m_D; 
% %V = diag(-1./x); 
% %V = V + diag(1i * (C_f*alpha_s*x.^2*Temp*m_D^2)/6 .*(2*gamma_E+(log(x.*m_D)).^2 - 8/3));
% 
% % Kinetic
% main_diagonal = (-2/(dx^2*M)) * ones(N, 1);
% off_diagonal = 1/(M*dx^2)*ones(N - 1, 1);
% T = diag(main_diagonal) + diag(off_diagonal, 1) + diag(off_diagonal, -1);
% 
% %Calculation
% H = -0.5*T + V;
% [Evec,Eval] = eig(H);
% Eval=diag(Eval);sort(-Eval);



%% Matrices and sparse (3D)Calculation
y = x;
z = x;
[xx,yy,zz] = meshgrid(x,y,z);
xx=xx(:);yy=yy(:);zz=zz(:);
I = eye(N); I = sparse(I);
%Potentials

%V = 0.5*diag(sqrt(xx.^2+yy.^2+zz.^2)); %Spherical Oscillator
V = -diag(1/sqrt(xx.^2+yy.^2+zz.^2)); % Coulomb

% Kinetic

main_diagonal = (-2/(dx^2*M)) * ones(N, 1);
off_diagonal = 1/(M*dx^2)*ones(N - 1, 1);
T = diag(main_diagonal) + diag(off_diagonal, 1) + diag(off_diagonal, -1);

V=sparse(V);
T=sparse(T);
Tx=kron(I,kron(I,T));
Ty=kron(I,kron(T,I));
Tz=kron(kron(T,I),I);
T = Tx + Ty + Tz;
H = -0.5*T+V;
clearvars T V Tx Ty Tz main_diagonal off_diagonal x y z I %delete vars to free memory

[Evec,Eval] = eigs(H,6,'sr');
Eval=diag(Eval);sort(-Eval);

%% Plots
%1D

%Real Potential

%Complex Potential

% subplot(2,2,1);
% plot(x,real(-Evec(:,1)),'.','markersize',12),grid on
% title(sprintf('Energy: %d = %20.13f',1,real(Eval(1))));
% 
% subplot(2,2,2);
% plot(x,imag(-Evec(:,1)),'.','markersize',12),grid on
% title(sprintf('Energy: %d = %20.13f',1,imag(Eval(1))));
% subplot(2,2,3);
% plot(x,real(-Evec(:,2)),'.','markersize',12),grid on
% title(sprintf('Energy: %d = %20.13f',2,real(Eval(2))));
% subplot(2,2,4);
% plot(x,imag(-Evec(:,2)),'.','markersize',12),grid on
% title(sprintf('Energy: %d = %20.13f',2,imag(Eval(2))));

%3D
EV1 = reshape((Evec(:,1)),N,N,N);
scatter3(xx(:),zz(:),yy(:),1,-EV1(:),'fill')