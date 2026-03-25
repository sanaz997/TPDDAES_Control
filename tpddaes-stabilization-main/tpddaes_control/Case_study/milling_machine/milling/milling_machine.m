function [A,B,C,T,E,tau_A,points_d]=milling_machine()
q=3/4;fz=0.1*1e-3;
phi_en=0;  % radiants
phi_ex=41.41*pi/180;      % radiants
ap=1.4*1e-3;
m=1;    omega=3200;% Kg
c=40;      % N*s/m
k=10^6;   % N/m
Kt=800*10^6;  % N/m^2
Kr=300*10^6;  % N/m^2
Z=2;  
T=60/Z/omega;
pp = ap*q*fz^(q-1);
qcoscos2=@(t) qcoscos(t,omega,Z,pp,q,phi_en,phi_ex);
qsincos2=@(t) qsincos(t,omega,Z,pp,q,phi_en,phi_ex);
points_d=disc_points(phi_en,phi_ex,omega,Z);
ka=23.3588;
%% 
% % %% In this part of the code, we have only considered the additioanl delay in the control input: 
E=[1 0  ;  c m  ];


A={@(t) [ 0 1 ; -k+Kt*qsincos2(t)-Kr*qcoscos2(t) 0  ];
   @(t) [0 0 ;-Kt*qsincos2(t)+Kr*qcoscos2(t) 0 ];};

B=[0 ;ka ];
C=[1 0];
tau_A={@(t)0 ; @(t) 60/Z/omega;};
