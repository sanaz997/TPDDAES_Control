function rddae=tpddae_create(A,B,C,E,tau_A,points_d,T)

% This function returns the struct rddae from the given Time-periodic delay
% differential equation of retarded type. 
% inputs:
% A: is defined as a cell-array where each cell is a time-dependent
% function handle d\times d. 
% B: constant matrix of size d\times r
% C: constant matrix of size p\times d
% E: constant matrix of size d\times d which is allowed to be singular
% T: time-period 
% tau_A: delay of the system in the form of cell array where each cell is
% allowed to be a time-dependent function handle.
% points_d: discontinuty times of the system in the interval [0,T] in the
% form of an array

rddae.A=A;
rddae.B=B;
rddae.E=E;
rddae.C=C;
rddae.T=T;
rddae.tau_A=tau_A;
rddae.points_d=points_d; 