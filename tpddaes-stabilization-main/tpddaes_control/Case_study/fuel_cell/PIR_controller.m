clc;clear all;close all;
%% add_path 
%%% this adds the folder with all its subfolders:
%%% Don't forget to change to your own path.
addpath(genpath('C:\Users\sanaz\OneDrive\Desktop\TPDDAEs_MATLAB'))
%%
%%% compare eigenvalue approximation via the two methods: 
[A,B,C,T,E,tau_A,points_d]=fuel_cell();
rddae=tpddae_create(A,B,C,E,tau_A,points_d,T);
kp=0.1;ki=1;kr=0.07;
gain.K=[-kp -ki kr];
gain.type='TC';
M=25;
[z,g_z,gg]=FM_computation(rddae,M,'gain',gain);
abs(z) %%% eigenvalue approximation using Pseduospectral collocation methods
cr=-16.3831; %%%% using TDS_Control software
z2=abs(exp(cr(1)*T)) %%% eigenvalue approximation using Floquet exponents
 %%
%  clear gain;
%  Theta=K_QR6.theta;
%  gain=K_QR6;
%  gain.K=[];
%  gain_tem.K=[];
% gain_mean.type='TC'
%  for i = 1:numel(fieldnames(Theta))
%      fieldName = sprintf('K_%d', i);
%     gain.K =[ gain.K  K_QR6.K.(fieldName)]
%  gain_tem.K=[ gain_tem.K ; K_QR6.K.(fieldName)]
%  end
%  gain_mean=gain;
%  gain_mean.K=[];
% gain_mean.K=mean(gain_tem.K)

%% 
[K_FR,z_FR,time_FR,pars]=tpddae_stabopt_static(rddae,M,'FR','n_FR',3,'nstart',3);
%%
close all;
[K_QR6,z_QR6,time_QR6,pars]=tpddae_stabopt_static(rddae,M,'QR','nstart',3,'lambda1',2e-5);
disp('Spectral radius with time-periodic gain')
abs(z_QR6)
disp('Spectral radius with time-constant gain')
abs(z2)
%%
close all;
figure;
subplot(1,3,1)
kp=0.1;ki=1;kr=0.07;
gain.K=[-kp -ki kr];
gain.type='TC';
load('PIR_FR.mat')
gain.M=M;
gain.T=T;
plot_K(gain)
legendInfo = {'\boldmath $\mathcal{K}_p$','\boldmath $\mathcal{K}_i$','\boldmath $\mathcal{K}_r$'};
lgd=legend(legendInfo, 'interpreter', 'latex');
ylim([-2 0.5])
title('TI')


subplot(1,3,2);
load('PIR_QR.mat')
K_QR6.points_d=T;
K_QR6.Delta=T;
plot_K(K_QR6);
legendInfo = {'\boldmath $\mathcal{K}_p$','\boldmath $\mathcal{K}_i$','\boldmath $\mathcal{K}_r$'};
lgd=legend(legendInfo, 'interpreter', 'latex');
ylim([-20 40])
title('QR')


subplot(1,3,3);
load('PIR_FR.mat')
plot_K(K_FR);
title('FR')
legendInfo = {'\boldmath $\mathcal{K}_p$','\boldmath $\mathcal{K}_i$','\boldmath $\mathcal{K}_r$'};
lgd=legend(legendInfo, 'interpreter', 'latex');
hold on;
ylim([-3 11])