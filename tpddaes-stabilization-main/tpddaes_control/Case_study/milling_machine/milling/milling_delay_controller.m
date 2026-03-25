clc;clear all;close all;

%%% test milling: 
M=40;
[A,B,C,T,E,tau_A,points_d]=milling_machine();
rddae=tpddae_create(A,B,C,E,tau_A,points_d,T);
[z,~]=FM_computation(rddae,M);

[K_FR,z_FR,time_FR,pars]=tpddae_stabopt_static(rddae,M,'FR','n_FR',5,'nstart',10,'normtol',1e-5);

z_FR_T0=z_FR;
K_FR_T0=K_FR;
pars_T0=pars;
time_FR_T0=time_FR;

%%
tau_out=T/300;
[A,B,C,T,E,tau_A,points_d]=milling_machine2(tau_out);
rddae=tpddae_create(A,B,C,E,tau_A,points_d,T);
[K_FR,z_FR,time_FR,pars]=tpddae_stabopt_static(rddae,M,'FR','n_FR',5,'nstart',10,'normtol',1e-5);
% 
% z_FR_T100=z_FR;
% K_FR_T100=K_FR;
% pars_T100=pars;
% time_FR_T100=time_FR;
%%
tau_out=T/100;
[A,B,C,T,E,tau_A,points_d]=milling_machine2(tau_out);
rddae=tpddae_create(A,B,C,E,tau_A,points_d,T);
[K_FR,z_FR,time_FR,pars]=tpddae_stabopt_static(rddae,M,'FR','n_FR',5,'nstart',10);

%%
tau_out=T/10;

[A,B,C,T,E,tau_A,points_d]=milling_machine2(tau_out);
rddae=tpddae_create(A,B,C,E,tau_A,points_d,T);
[K_FR,z_FR,time_FR,pars]=tpddae_stabopt_static(rddae,M,'FR','n_FR',5,'nstart',10);
%%
close all;
figure();
addpath(genpath('C::\Users\sanaz\OneDrive\Desktop\TPDDAEs_MATLAB'))
%%
close all;
load('T0.mat')
hold on
subplot(2,2,1);
K_FR.Delta=T;
K_FR.points_d=points_d;

plot_K(K_FR)
xlim([0 T])
legend('\boldmath$\mathcal{K}_1$', 'interpreter', 'latex')
xlabel('t(s)')
title('\boldmath $ \tau_\zeta=0$', 'interpreter', 'latex')
hold on;
%%
hold on
subplot(2,2,2);
load('T300.mat')
K_FR_T300.Delta=T;
K_FR_T300.points_d=points_d;
plot_K(K_FR_T300)
xlim([0 T])
legend('\boldmath$\mathcal{K}_1$','\boldmath$\mathcal{K}_2$', 'interpreter', 'latex')
xlabel('t(s)')
title('\boldmath $ \tau_\zeta={T}/{300}$', 'interpreter', 'latex')
hold on;
%%
hold on
subplot(2,2,3);
load('T100.mat')
K_FR_T100.Delta=T;
K_FR_T100.points_d=points_d;
plot_K(K_FR_T100)
xlim([0 T])
legend('\boldmath$\mathcal{K}_1$','\boldmath$\mathcal{K}_2$', 'interpreter', 'latex')
xlabel('t(s)')
title('\boldmath $ \tau_\zeta={T}/{100}$', 'interpreter', 'latex')
hold on;
%%
hold on
subplot(2,2,4);
load('T10.mat')
K_FR_T10.Delta=T;
K_FR_T10.points_d=points_d;
plot_K(K_FR_T10)
xlim([0 T])
legend('\boldmath$\mathcal{K}_1$','\boldmath$\mathcal{K}_2$', 'interpreter', 'latex')
xlabel('t(s)')
title('\boldmath $ \tau_\zeta={T}/{10}$', 'interpreter', 'latex')
hold on;
%%
subplot(2,2,1)
ylim([-8e4,0e4])
subplot(2,2,2)
ylim([-1.5e5 ,1.5e5])
subplot(2,2,3)
ylim([-1.5e5 ,1.5e5])
subplot(2,2,4)
ylim([-6e4 ,6e4])





