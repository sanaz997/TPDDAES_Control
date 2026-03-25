clc;clear all;close all;
%%
%%%%%
% This file studies the effect of phase mismatch of the time-periodic feedback gain and the time-periodic
% system. this case study is for delayed Mathieu equation.
%% add_path 
%%% this adds the folder with all its subfolders:
addpath(genpath('C:\Users\sanaz\OneDrive\Desktop\TPDDAEs_MATLAB'))
[A,B,C,T,E,tau_A,points_d]=delayed_mathieu();
M=28;
rddae=tpddae_create(A,B,C,E,tau_A,points_d,T);
[z,~]=FM_computation(rddae,M);
abs(z);
%%
clear gain
load('mathieu_stat.mat');
gain.K=[K_FR5.K_1 K_FR5.K_2 K_FR5.K_3];
gain.type='FR';
gain.n_FR=3;
[z,~]=FM_computation(rddae,M,'gain',gain);
abs(z);
%%
clear h_vec z3_vec opt_z
close all
AA=@(t) rddae.A{1}(t);
h=0;
hmax=2*pi;
end_hmax=2*hmax;
i=1;
for h=0:end_hmax/200:end_hmax
A{1}=@(t) AA(t)+B*(K_FR5.K_1+K_FR5.K_2.*cos(2*pi*t/T-h)+K_FR5.K_3.*sin(2*pi*t/T-h))*C;
rddae2=tpddae_create(A,B,C,E,tau_A,points_d,T);
[z3,g]=FM_computation(rddae2,M);
abs(z3);
h_vec(i,1)=h;
z3_vec(i,1)=abs(z3);
opt_z(i,1)=abs(z);
z_TC(i,1)=1.1824;
z_stab(i,1)=1;
i=1+i;
end
%%
close all;
figure();
plot(h_vec,z_TC,':','LineWidth',5,'Color', [0.9290 0.6940 0.1250])
hold on
plot(h_vec,opt_z,'-.','LineWidth',5,'Color', [0.8500 0.3250 0.0980])


hold on

hold on
plot(h_vec, z3_vec,'LineWidth',5,'Color', [0 0.4470 0.7410])
yline(1, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 4);
xlim([0 4*pi])
xticks([0, pi, 2*pi, 3*pi, 4*pi]);
grid on;
% Set x-axis tick labels as multiples of pi
xticklabels({'0', '\pi', '2\pi', '3\pi','4\pi'});
ylabel('\boldmath $\mathcal{P}(\mathcal{U}_M)$','interpreter', 'latex' ,'FontSize',16, 'FontWeight','bold' )
xlabel('\boldmath $\Phi (\mathrm{rad})$' ,'FontSize',16,'interpreter', 'latex', 'FontWeight','bold')
legend( '\boldmath $\mathcal{P}\big(\mathcal{U}_M(K)\big)$','\boldmath $\mathcal{P}\big(\mathcal{U}_M(\mathcal{K}(t))\big)$ \textbf{with} \boldmath $\mathcal{K}(t)=\mathcal{K}(t+T)$ ','\boldmath $\mathcal{P}(\mathcal{U}_M(\mathcal{K}(t-\Phi))\big)$ \textbf{with} \boldmath $\mathcal{K}(t)=\mathcal{K}(t+T)$','FontSize',25, 'interpreter', 'latex')

ax = gca;
ax.FontSize =25;  
ax.FontWeight = 'bold';
ax.LineWidth = 2;
grid on
ax.FontWeight = 'bold';
set(ax, 'LineWidth', 3);
ax = gca;
ax.FontSize = 25;  
set(gca, 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

