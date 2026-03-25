clc;clear all;close all;
%% add_path 
%%% this adds the folder with all its subfolders:
addpath(genpath('C:\Users\sanaz\OneDrive\Desktop\TPDDAEs_MATLAB'))
%%
%%%% Time-periodic linearized dynamics:
%%%% system parameters:
lambda=-0.02;gamma=1;
omega=1-lambda*gamma;T=2*pi/omega;
%%%%% periodic solution:
x_1=@(t) sqrt(-lambda)*cos(omega*t);
x_2=@(t) sqrt(-lambda)*sin(omega*t);
A_0=@(t) [lambda+3*x_1(t)^2+x_2(t)^2-2*gamma*x_1(t)*x_2(t)   -1+2*x_1(t)*x_2(t)-gamma*x_1(t)^2-3*gamma*x_2(t)^2;
  1+3*gamma*x_1(t)^2+gamma*x_2(t)^2+2*x_1(t)*x_2(t)   lambda+2*gamma*x_1(t)*x_2(t)+x_1(t)^2+3*x_2(t)^2  
];
b_vec=[0;1];
%%
%%%%%%%% TPDDAEs Formulation:
E=[eye(2,2) zeros(2,2);zeros(2,4)];
A={@(t) [A_0(t)  zeros(2,2);zeros(2,2) eye(2,2)] ;@(t) [zeros(2,4); -eye(2,2) zeros(2,2) ]};
B=[b_vec;zeros(2,1)];
C=[-eye(2,2) eye(2,2)];
tau_A={@(t) 0; @(t) T};
points_d=[];
%%
figure();
 rddae=tpddae_create(A,B,C,E,tau_A,points_d,T);
 M=45;
 clc
 gain.type='TC'
 gain.K=[0 0];
[z,g_z,gg,d_vec]=FM_computation(rddae,M);


figure;
plot(real(d_vec), imag(d_vec), 'ro', 'MarkerSize', 8, 'DisplayName', 'Eigenvalues');
hold on;

% Plot the unit circle
theta = linspace(0, 2*pi, 1000);
unit_circle = exp(1i * theta); % Points on the unit circle
plot(real(unit_circle), imag(unit_circle), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Unit Circle');

% Configure plot
axis equal;
grid on;
xlabel('Real Part');
ylabel('Imaginary Part');
title('Eigenvalues and Unit Circle');
legend('show');
hold off;
%%%%%%%%%%
exp(-2*lambda*T)
abs(z)

%%
load('K_FR22.mat')
K_FR.K=[ K_FR.K_1 K_FR.K_2 K_FR.K_3 K_FR.K_4 K_FR.K_5];
K_FR.T=T;
K_FR.n_FR=5;
K_FR.type='FR';
[z,g_z,gg,d_vec2]=FM_computation(rddae,M,'gain',K_FR);

% Plot the eigenvalues
hold on 
plot(real(d_vec2), imag(d_vec2), 'ro', 'color','g','MarkerSize', 8, 'DisplayName', 'Eigenvalues');
hold on;
d_vec3=sort(abs(d_vec2));
d_vec3(end-5:end,1)
%%
M=45;
[K_FR,z_FRR,time_FR,pars]=tpddae_stabopt_UPO(rddae,M,'FR','n_FR',5,'nstart',3);
%%
figure();
plot_K(K_FR)
%%
%%%%%% Here we check the behavior of the nonlinear system with the Pyragas
%%%%%% type-feedback controller:
clc;close all
lags=[T];%%% time-delay
 tspan=[0,8*T]; %%% time-interval
options = ddeset('RelTol',1e-4,'AbsTol',1e-4, 'MaxStep', 1e-3);
sol_ideal = dde23(@(t, x, x_del) TPDDEs_op_nl(t, x, x_del),lags,@(t) initial(t,0),tspan,options);
 epsilon=5e-3;
sol_open_per = dde23(@(t, x, x_del) TPDDEs_op_nl(t, x, x_del),lags,@(t) initial(t,epsilon),tspan,options);
 hold on;
sol_clos_init = dde23(@(t, x, x_del) TPDDEs_nl(t, x, x_del, K_FR),lags,@(t) initial(t,epsilon),tspan,options);
lambda=-0.02;
close all
clear norm_vec sol ideal_vec sol2 norm_vec2
sol=sol_clos_init;
sol2=sol_open_per;
%%
figure;
for i=1:length(sol.y)
norm_vec(i,1)=abs(norm(sol.y(:,i)));
ideal_vec(i,1)=sqrt(-lambda);
end
plot(sol.x,ideal_vec,'--','LineWidth',5)
hold on;
plot(sol.x,norm_vec,'LineWidth',5)
hold on; 


for i=1:length(sol2.y)
norm_vec2(i,1)=abs(norm(sol2.y(:,i)));
end
plot(sol2.x,norm_vec2,':','LineWidth',5)
grid on;
legendInfo=legend('UPO','With Controller','Without Controller');
ylim([0.05,0.25])
xlim([0 8*T])


grid on;
xlabel('t(s)','FontSize',25, 'FontWeight','bold');
ylabel(' \boldmath $\|x(t)\|$ ','FontSize',25, 'interpreter', 'latex');
ax = gca;
ax.FontSize = 20;  
set(ax, 'LineWidth', 3);
ax = gca;
ax.FontSize = 30;  
set(gca, 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

%%
figure();
plot_K(K_FR)
%%
figure;
plot(sol_ideal.y(1,:),sol_ideal.y(2,:),'--','LineWidth',2);
hold on;
plot(sol_open_per.y(1,:),sol_open_per.y(2,:),'LineWidth',2);
hold on;
plot(sol_clos_init.y(1,:),sol_clos_init.y(2,:),'LineWidth',2);
legend('open_loop on the periodic orbit','open-loop in the neighbour' ,'closed-loop in the neighbor')

plot(sol_clos_init.y(1,1),sol_clos_init.y(2,1), 'o', 'MarkerFaceColor', [0.7 0 0], 'MarkerEdgeColor', [0.7 0 0], 'MarkerSize', 12, 'Color','#A2142F')
text(sol_clos_init.y(1,1),sol_clos_init.y(2,1), sprintf('initial'), 'FontSize', 25, 'FontWeight', 'bold');
grid on;
hold on;
%%
function dxdt =TPDDEs_op_nl(t,x,x_del)
lambda=-0.02;gamma=1;omega=1-lambda*gamma;T=2*pi/omega;
%%%%% periodic solution:
F=@(t)[
  lambda*x(1)-x(2)+(x(1)^2+x(2)^2)*(x(1)-gamma*x(2));
  x(1)+lambda*x(2)+(x(1)^2+x(2)^2)*(gamma*x(1)+x(2))
];
dxdt=F(t);
end
%%
function dxdt =TPDDEs_nl(t,x,x_del,K_FR)
lambda=-0.02;gamma=1;omega=1-lambda*gamma;T=2*pi/omega;
%%%% periodic solution:
F=@(t)[
  lambda*x(1)-x(2)+(x(1)^2+x(2)^2)*(x(1)-gamma*x(2));
  x(1)+lambda*x(2)+(x(1)^2+x(2)^2)*(gamma*x(1)+x(2))
];
b_vec=[0;1];
K_2=@(t) K_FR.K_1*1 +K_FR.K_2*cos(omega*t)+ K_FR.K_3*sin(omega*t)+K_FR.K_4*cos(omega*t/2)+ K_FR.K_5*sin(omega*t/2);
Pyr=@(t) b_vec*K_2(t)*(x_del-x);
 dxdt=F(t)+Pyr(t);
end
%%
function out=initial(t,epsilon)
lambda=-0.02;gamma=1;
omega=1-lambda*gamma;T=2*pi/omega;
d_x_1=@(t) sqrt(-lambda)*cos(omega*t)+epsilon; %%%%% assuming that the paper of smapled data approach made mistake :-
d_x_2=@(t) sqrt(-lambda)*sin(omega*t)+1/8*epsilon;
out=[d_x_1(t);d_x_2(t)];
end




