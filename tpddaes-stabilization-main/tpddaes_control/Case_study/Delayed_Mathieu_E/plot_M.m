 clc;clear all;close all;
%%%%% This file plots the relative error of the spectral radius approximation with respect to 
% different choices of the discretization degree%% add_path 
%%
%%% This adds the folder along with all its subfolders.
addpath(genpath('C:\Users\sanaz\OneDrive\Desktop\TPDDAEs_MATLAB'))
 %% 
[A,B,C,T,E,tau_A,points_d]=delayed_mathieu();
rddae=tpddae_create(A,B,C,E,tau_A,points_d,T);
[A,B,C,T,E,tau_A,points_d]=delayed_mathieu();
M_final=70;
M_init=1;
j=1;
[z_final,g_z]=FM_computation(rddae,M_final);

for M=M_init:1:M_final
    [z,g_z]=FM_computation(rddae,M);
    z_vec(1,j)=z;  
    M_vec(1,j)=M;
    err(1,j)=norm((abs(z_vec(1,j))-abs(z_final))/abs(z_final));
    j=j+1;
end

    %%
    figure();
%%
%%%%%%%%%% You can comment this out if you do not want to load the data from the workspace
load('fig_M.mat')
%%
    semilogy(M_vec(1:end-1),err(1:end-1),'LineWidth',5);
    grid off;
    grid on;
    
    ylabel('Relative Error','FontSize',25,'FontWeight','bold')
    xlabel('Dsicretization Degree ','FontSize',25,'FontWeight','bold')
    title('Relative error in approximating the spectral radius with respect to M','FontSize',16,'FontWeight','bold')
    ax = gca;
    ax.FontSize =25;  
    xlim([M_init M_final-1])
    
    
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
    

