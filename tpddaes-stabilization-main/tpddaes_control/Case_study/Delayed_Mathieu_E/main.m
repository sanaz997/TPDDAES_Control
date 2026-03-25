 clc;clear all;close all;
%% add_path 
%%% this adds the folder with all its subfolders:
addpath(genpath('C:\Users\sanaz\OneDrive\Desktop\TPDDAEs_MATLAB'))

 %%
%%%% delayed Mathieu equations: 
[A,B,C,T,E,tau_A,points_d]=delayed_mathieu();
M=28;
rddae=tpddae_create(A,B,C,E,tau_A,points_d,T);
[z,~]=FM_computation(rddae,M);
abs(z)
%%%%% Static Feedback Controller with Time_constant Feedback Gain
%%
clc
[K_TC,z_TC2,time_TC]=tpddae_stabopt_static(rddae,M,'TC','nstart',3);
[z2,~]=FM_computation(rddae,M,'gain',K_TC);
abs(z2)
plot_K(K_TC)
 %%
%%%%% Static Feedback Controller with Time_periodic Feedback Gain:
%% QR approach:
[K_QR6,z_QR6,time_QR6]=tpddae_stabopt_static(rddae,M,'QR','nstart',3,'lambda1',7e-4);
%%
plot_K(K_QR6)
%% FR approach:
[K_FR5,z_FR5,time_FR5]=tpddae_stabopt_static(rddae,M,'FR','nstart',3,'n_FR',3);
plot_K(K_FR5)
%%
%%%%% plot_figures of the paper with the data from the workspace
clc;clear all;close all
load('mathieu_stat.mat')
figure();
 subplot(1,3,1)
K_TC.points_d=points_d;
K_TC.Delta=T;
plot_K(K_TC);
ylim([-6 -2])
hold on; 
subplot(1,3,2)
K_QR.points_d=points_d;
K_QR.Delta=T;
plot_K(K_QR)
hold on;
subplot(1,3,3)
plot_K(K_FR5);
%%% 
%% Dynamic Controller Design:
%%%%% Dynamic Feedback Controller with Time_constant Feedback Gain:
q=2;
[rddae_2,K_TC_dyn,z_TC_dyn,time_TC_dyn,pars,gainT_dyn]=tpddae_stabopt_dynamic(rddae,q,M,'TC','nstart',3);
%%
gainT_dyn.M=M;
plot_K_dynamic(rddae_2,q,gainT_dyn,pars.slack)
 %%
%%%%% plot_figures of the paper with the data from the workspace
clc;clear all;close all;
load('dyn_TC.mat')
plot_K_dynamic(rddae_2,q,gainT_dyn,pars.slack)
subplot(2,2,1);
ylim([-5,13])
subplot(2,2,2)
ylim([-4, 2])
subplot(2,2,3)
ylim([-4 ,4])
subplot(2,2,4)
ylim([-5,8])
 %%
 %%% QR apprach:
 q=2;
 [rddae_2,K_QR_dyn,z_QR_dyn,time_QR_dyn,pars,gainQ_dyn]=tpddae_stabopt_dynamic(rddae,q,M,'QR','nstart',3,'lambda1',2e-5);
%%
clc;clear all;close all;
%%%%% plot_figures of the paper with the data from the workspace
load('dyn_QR.mat')
% gainQ_dyn.points_d=points_d;
% gainQ_dyn.Delta=T;
plot_K_dynamic(rddae_2,q,gainQ_dyn,pars.slack,pars.slack_mat)
subplot(2,2,1);
ylim([-5,5])
subplot(2,2,2)
ylim([-4,5])
subplot(2,2,3)
ylim([-3,3])
subplot(2,2,4)
ylim([-15,5])
%%
q=2;
[rddae_2,K_FR_dyn,z_FR_dyn,time_FR_dyn,pars,gainF_dyn]=tpddae_stabopt_dynamic(rddae,q,M,'FR','nstart',10,'n_FR',3);
%%
clc;clear all;close all;
%%%%% plot_figures of the paper with the data from the workspace
load('work_space\dyn_FR.mat');
plot_K_dynamic(rddae_2,q,gainF_dyn,pars.slack);
subplot(2,2,1);
ylim([-10,7])
subplot(2,2,2)
ylim([-6,6])
subplot(2,2,3)
ylim([-7,5])
subplot(2,2,4)
ylim([-10,5])