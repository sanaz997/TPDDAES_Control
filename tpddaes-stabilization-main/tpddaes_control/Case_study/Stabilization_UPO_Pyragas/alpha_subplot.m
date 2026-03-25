clc;clear all;close all;
number=4;
load('K_FR22.mat')
K=[ K_FR.K_1 K_FR.K_2 K_FR.K_3 K_FR.K_4 K_FR.K_5];
[z,g_z,gg,d_vec]=FM_computation(rddae,M);
modulus = abs(d_vec);
[~, sorted_indices] = sort(modulus, 'descend');
top_5_indices = sorted_indices(1:number);
d_vec_top_init = d_vec(top_5_indices);
figure;
theta = linspace(0, 2*pi, 1000);
unit_circle = exp(1i * theta); 
h1=plot(real(unit_circle), imag(unit_circle), 'b', 'LineWidth', 3);
hold on;
eig_plot(d_vec_top_init,'r')
hold on;
inv_zero= 1.000000000000001e+00 + 0.000000000000000e+00i;
h2=plot(real(inv_zero), imag(inv_zero), "p",'Color','r', 'MarkerSize', 25, 'LineWidth', 3);
unstable_eig= 1.279410498514245e+00 + 0.000000000000000e+00i;
h3=plot(real(unstable_eig), imag(unstable_eig), "o",'Color','r', 'MarkerSize', 25, 'LineWidth', 3);
K_FR.K=K;
[z,g_z,gg,d_vec2]=FM_computation(rddae,M,'gain',K_FR);
modulus = abs(d_vec2);

stable_eig=1.210776278161328e-01 + 0.000000000000000e+00i;
h4=plot(real(stable_eig), imag(stable_eig), "X",'Color',[0.4660 0.6740 0.1880], 'MarkerSize', 25, 'LineWidth', 3);
[~, sorted_indices] = sort(modulus, 'descend');
top_5_indices = sorted_indices(1:number);


%%
grid on;
d_vec2_top= d_vec2(top_5_indices);
eig_plot_cl(d_vec2_top,'r')
%%
%%
i=1;
for alpha=0:1e-1:1
    K_FR.K=alpha*K;
[z,g_z,gg,d_vec2]=FM_computation(rddae,M,'gain',K_FR);
modulus = abs(d_vec2);
[~, sorted_indices] = sort(modulus, 'descend');
top_5_indices = sorted_indices(1:number);
d_vec2_top= d_vec2(top_5_indices);
threshold = 1e-9;
modulus = abs(d_vec2_top);
is_almost_one = abs(modulus - 1) < threshold; 
filtered_eigenvalues = d_vec2_top(~is_almost_one);
filtered_eigenvalues;
Map(:,i)=filtered_eigenvalues;
i=i+1;
end
%%
%%
hold on;
eig_plot_map(Map)
hold on; 
xlim([-0.2 1.4])
ylim([-0.3 0.3])
h2=plot(real(inv_zero), imag(inv_zero), "p",'Color','r', 'MarkerSize', 25, 'LineWidth', 3);
legend([h1,h2,h3,h4],'Unit Disk','Invariant Eigenvalue','Open-Loop Eigenvalues','Closed-Loop Eigenvalues', 'latex', 'FontWeight', 'bold')
ax = gca;
ax.FontWeight = 'bold';
ax.LineWidth = 3;
grid on
ylabel('\boldmath{$\mathrm{Im}(\mu)$}','FontSize',25, 'interpreter', 'latex', 'FontWeight', 'bold')
xlabel('\boldmath{$\mathrm{Re}(\mu)$}','FontSize',25, 'interpreter', 'latex', 'FontWeight', 'bold')
ax.FontWeight = 'bold';
set(ax, 'LineWidth', 3);
ax = gca;
ax.FontSize = 30;  
set(ax, 'LineWidth', 3);
ax = gca;
ax.FontSize = 30;  
set(gca, 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');



%%
function eig_plot(d_vec,color)
tolerance=1e-9;
for i=1:length(d_vec)
    d_vec(i)
if abs(d_vec(i))>1+tolerance
hold on;
elseif abs(d_vec(i))<1-tolerance
plot(real(d_vec(i)), imag(d_vec(i)), 'ro','Color','r', 'MarkerSize', 25, 'LineWidth', 3);
hold on;
end
end
end





%%
function eig_plot_cl(d_vec,color)
tolerance=1e-9;
for i=1:length(d_vec)
    
if abs(d_vec(i))>1+tolerance
plot(real(d_vec(i)), imag(d_vec(i)), '^','Color',[0.4660 0.6740 0.1880], 'MarkerSize', 8);
hold on;
elseif abs(d_vec(i))<1-tolerance
 d_vec(i)
plot(real(d_vec(i)), imag(d_vec(i)), 'X','Color',[0.4660 0.6740 0.1880], 'MarkerSize', 26, 'LineWidth', 3);
hold on;
end
end
end
