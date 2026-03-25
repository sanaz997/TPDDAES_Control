function [A,B,C,T,E,tau_A,points_d]=milling_machine2(tau_out)
q=3/4;fz=0.1*1e-3;
phi_en=0;  % radiants
phi_ex=41.41*pi/180;      % radiants
ap=1.4*1e-3;
m=1;    omega=3200;% Kg
%%%omega=5000;
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

E=[1 0  0;  c m  0;   0 0 0  ];

A={@(t) [ 0 1 0;  -k+Kt*qsincos2(t)-Kr*qcoscos2(t) 0 0 ; 0  0 -1];
        @(t) [0 0 0; -Kt*qsincos2(t)+Kr*qcoscos2(t) 0 0;0 0 0  ];
 @(t) [0 0 0;    0 0 0; 1 0 0 ];};

B=[0;ka;0 ];
C=[1 0 0; 0 0 1 ];
tau_A={@(t)0 ; @(t)60/Z/omega;@(t) tau_out;};
