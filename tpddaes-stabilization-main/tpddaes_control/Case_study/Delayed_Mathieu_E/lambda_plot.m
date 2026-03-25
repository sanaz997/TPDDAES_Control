clc;clear all;close all;
%%
%%% This file studies the balance between the regularization term and the
% objective function for different values of the penalty parameter $\lambda_1$
%% add_path 
%%% this adds the folder with all its subfolders:
addpath(genpath('C:\Users\sanaz\OneDrive\Desktop\TPDDAEs_MATLAB'))
%%
num=6; %%% number of scenarios
lambda0=1e-2; %%% initial and highest value of lambda1
lambda1=lambda0;
[A,B,C,T,E,tau_A,points_d]=delayed_mathieu();
M=28;
rddae=tpddae_create(A,B,C,E,tau_A,points_d,T);
[z,g_z,gg]=FM_computation(rddae,M);
abs(z)
%%%%% Static Feedback Controller with Time_constant Feedback Gain:
clc
[K_opt,z_opt,elapsed_time,pars,obj,reg]=tpddae_stabopt_static(rddae,M,'QR','nstart',10,'lambda1',lambda0);

%%
i=1;
F11(i,1)=abs(z_opt); %%% Calculation of the objective function
F22(i,1)=abs(reg);%%% Caclulation of the newly defined regularization term
lambda_vecc(i,1)=lambda1;
%
update=1e-1;

for i=2:num

    lambda1=lambda1*update;
    [K_opt,z_opt,elapsed_time,pars,obj,reg]=tpddae_stabopt_static(rddae,M,'QR','nstart',10,'lambda1',lambda1);
    F11(i,1)=abs(z_opt) ;%%% Calculation of the objective function
    F22(i,1)=abs(reg);%%% Caclulation of the regularization term
    lambda_vecc(i,1)=lambda1;

end

%%
figure();
%%%%% If you want to obtain the figure of the paper uncomment the following
%%%%% load command:
load('mathieu_stat.mat')

plot(F11, F22,'LineWidth',5, 'Color','#0072BD');
hold on
plot(F11, F22, 'o', 'MarkerFaceColor', [0.7 0 0], 'MarkerEdgeColor', [0.7 0 0], 'MarkerSize', 12, 'Color','#A2142F'); % filled circles to denote the scanrios

for i = 1:length(F11)
    text(F11(i, 1)+2e-2,F22(i, 1)+5e2, sprintf('\\gamma_1=%.3e', lambda_vecc(i,1)), 'FontSize', 25, 'FontWeight', 'bold');
end

grid on;
xlabel('\boldmath $\mathcal{P}(\mathcal{U}_M) $','FontSize',25, 'interpreter', 'latex', 'FontWeight', 'bold')
ylabel('\boldmath $\mathrm{Tr}(K^T\Lambda_1K)$ ','FontSize',20, 'interpreter', 'latex', 'FontWeight', 'bold')
%title('Comparison the trade-off between the spectral radius and the regularization term','FontSize',25, 'FontWeight', 'bold')
ax = gca;

set(ax, 'LineWidth', 3);
ax = gca;
ax.FontSize = 30;  
set(gca, 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
xlim([0.1 1.4])


