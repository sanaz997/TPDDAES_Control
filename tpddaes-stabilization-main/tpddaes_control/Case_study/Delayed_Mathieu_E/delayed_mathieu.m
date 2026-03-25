function [A,B,C,T,E,tau_A,points_d]=delayed_mathieu()

delta=-4; epsilon=1;
omega=2*pi; 
E=[1 0 0;0 1 0;0 0 0];
A_0=@(t) [0 1 0 ; -(delta+epsilon*cos(omega*t)) 0 0 ; 0 0 -1];
A_1=@(t) [0 0 0 ;0 0 1 ;0 0 0];
B=[0;0;1];
C=[1 0 0;0 1 0];
A={@(t) A_0(t); @(t) A_1(t)};
T=2*pi/omega;  
tau_u=4*T/5;
tau_A={@(t) 0;@(t) tau_u};
points_d=[]; %%% discontinuity points
end