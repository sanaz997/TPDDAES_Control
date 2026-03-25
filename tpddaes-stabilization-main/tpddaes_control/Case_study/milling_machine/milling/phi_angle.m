function [ phi ] = phi_angle( j,t,omega,N )
phi=  2*pi*mod(omega*t/60+j/N,1);
end