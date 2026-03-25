clc;clear all;close all;

q=3/4;fz=0.1*1e-3;
phi_en=0;  % radiants
phi_ex=41.41*pi/180;      % radiants
ap=1.4*1e-3;  
omega=3200;% Kg
%%%omega=5000;
c=40;      % N*s/m
k=10^6;   % N/m
Z=2;  
T=60/Z/omega;
pp = ap*q*fz^(q-1);
qcoscos2=@(t) qcoscos(t,omega,Z,pp,q,phi_en,phi_ex);
qsincos2=@(t) qsincos(t,omega,Z,pp,q,phi_en,phi_ex);
N=2;Kt=800*10^6;  % N/m^2
Kr=300*10^6;  % N/m^2
m=1;  




i=1;
for t=0:1e-6:2*60/(N*omega)
   t_vec(i) =t;
   G_vec(i)=-Kt*qsincos2(t)+Kr*qcoscos2(t);
   i=1+i;
end

plot(t_vec,G_vec,'LineWidth',5, 'Color','#0072BD')







ax = gca;
ax.FontWeight = 'bold';
ax.LineWidth = 2;
grid on
xlabel('\textbf{t(s)}','FontSize',35, 'interpreter', 'latex', 'FontWeight', 'bold')
ylabel('\boldmath{$\tau$}-\textbf{periodic Function}\boldmath $~G (N/m)$','FontSize',35, 'interpreter', 'latex', 'FontWeight', 'bold')
ax.FontWeight = 'bold';
set(ax, 'LineWidth', 3);
ax = gca;
ax.FontSize = 30;  
set(gca, 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
xlim([0 2*T])

