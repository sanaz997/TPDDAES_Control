function plot_K_dynamic(rddae,q,K_fun,slack_out,varargin)
% This function plots the feedback gain value with respect to time over the interval [0,T]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=K_fun.M;
B=rddae.B;
C=rddae.C;
[r]=size(B,2);
[p]=size(C,1);
T=K_fun.T;
 intR=slack_out.intR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=100;
vec_x=zeros(T/(T/h)+1,1);
val=zeros(r*p,T/(T/h)+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if K_fun.type=='QR' %%% Feedback gain designed based on total variation method
    slack_mat=varargin{1};
%%
    x=1;
    for j=0:T/h:T    
        gain2=K_continous(K_fun,j,intR,slack_mat);
        val(1:r*p,x)=gain2(:); %%% Value of the feedback gain
        vec_x(x,1)=j; %%% Time vector
        x=x+1;
    end
    
    figure();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif K_fun.type=='TC' %%% Feedback gain designed based on truncated Fourier method

    gain=K_fun.K;
    
    x=1;
    for j=0:T/h:T    
    val(1:r*p,x)=gain(:); %%% Value of the feedback gain
    vec_x(x,1)=j; %%% Time vector
    x=x+1;
    end
    
    figure();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif K_fun.type=='FR' %%% Feedback gain designed based on truncated Fourier method
    n_FR=K_fun.n_FR;
    
    
    for i = 1:n_FR
         fieldName = sprintf('K_%d', i);
        gain_(1:r,p*(i-1)+1:p*i)=K_fun.(fieldName);  %%% storing feedback gain from the structure array
    end
    
    
    
    
    FR_matrix(1,1)=1;
     x=1;
     for j=0:T/h:T    
    
         for s=2:2:n_FR
            FR_matrix(s,1)=cos(2*pi*s/2*j/T);
            FR_matrix(s+1,1)=sin(2*pi*s/2*j/T);
         end
    
    FR_matrix2=kron(FR_matrix,eye(p));
    sum_val2=(gain_)*FR_matrix2;
    val(1:r*p,x)=sum_val2(:); %%% Value of the feedback gain
    vec_x(x,1)=j; %%% Time vector
    x=x+1;
    end
    
    figure();
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r=r-q;
    p =p-q;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,2,1)
    for k=1:q^2
        plot(vec_x,val(k,:),'LineWidth',5)
        hold on
    end

grid on;
xlabel('t(s)','FontSize',20,'FontWeight', 'bold')
ylabel(' \boldmath $\mathcal{A}_c(t)$ ','FontSize',20, 'interpreter', 'latex', 'FontWeight', 'bold')
%%% legend:
s=1;
for j = 1:q
    for i=1:q
     legendInfo{s} = ['$\boldmath \mathcal{A}_{\boldmath c_{' num2str(i), num2str(j),'}}$']; 
    s=s+1;
    end
end

lgd=legend(legendInfo, 'interpreter', 'latex');
fontsize(lgd,20,'points');
ax = gca;
ax.FontSize = 20;  
xlim([0 T])
ax.FontWeight = 'bold';
ax.LineWidth = 2;
grid on
ax.FontWeight = 'bold';
set(ax, 'LineWidth', 3);
ax = gca;
ax.FontSize = 20;  
set(gca, 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear legendInfo
    hold on
    
    subplot(2,2,2)
     for k=q^2+1:q^2+q*p
         plot(vec_x,val(k,:),'LineWidth',5)
         hold on
     end


grid on;
xlabel(' t(s)','FontSize',20,'FontWeight', 'bold')
ylabel(' \boldmath $\mathcal{B}_c(t)$ ','FontSize',20, 'interpreter', 'latex', 'FontWeight', 'bold')
%%% legend:
s=1;
for j = 1:p
    for i=1:q
    legendInfo{s} = ['$\boldmath \mathcal{B}_{\boldmath c_{' num2str(i), num2str(j),'}}$']; 
    s=s+1;
    end
end

lgd=legend(legendInfo,'interpreter', 'latex');
fontsize(lgd,20,'points');
ax = gca;
ax.FontSize = 20;  
xlim([0 T])
ax.FontWeight = 'bold';
ax.LineWidth = 2;
grid on
ax.FontWeight = 'bold';
set(ax, 'LineWidth', 3);
ax = gca;
ax.FontSize = 20;  
set(gca, 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear legendInfo
hold on

subplot(2,2,3)
 for k=q^2+q*p+1:q^2+q*p+q*r
     plot(vec_x,val(k,:),'LineWidth',5)
     hold on
 end

grid on;
xlabel('t(s)','FontSize',20,'FontWeight', 'bold')
ylabel(' \boldmath $\mathcal{C}_c(t)$ ','FontSize',20, 'interpreter', 'latex', 'FontWeight', 'bold')

%%% legend:
s=1;
for j = 1:q
    for i=1:r
    legendInfo{s} = ['$\boldmath \mathcal{C}_{\boldmath c_{' num2str(i), num2str(j),'}}$']; 
    s=s+1;
    end
end

lgd=legend(legendInfo, 'interpreter', 'latex');
fontsize(lgd,20,'points');
ax = gca;
ax.FontSize = 20;  
xlim([0 T])
ax.FontWeight = 'bold';
ax.LineWidth = 2;
grid on
ax.FontWeight = 'bold';
set(ax, 'LineWidth', 3);
ax = gca;
ax.FontSize = 20;  
set(gca, 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear legendInfo
hold on
%ylim([-6 4])
subplot(2,2,4)
 for k=q^2+q*p+q*r+1:q^2+q*p+q*r+r*p
     plot(vec_x,val(k,:),'LineWidth',5)
     hold on
 end

 
 grid on;
xlabel('t(s) ','FontSize',20,'FontWeight', 'bold')
ylabel(' \boldmath $\mathcal{D}_c(t)$ ','FontSize',20, 'interpreter', 'latex', 'FontWeight', 'bold')
%%% legend:
s=1;
for j = 1:p
    for i=1:r
      legendInfo{s} = ['$\boldmath \mathcal{D}_{\boldmath c_{' num2str(i), num2str(j),'}}$']; 
    s=s+1;
    end
end

lgd = legend(legendInfo, 'Interpreter', 'latex');
fontsize(lgd,20,'points');
ax = gca;
ax.FontSize = 20;  
xlim([0 T])
ax.FontWeight = 'bold';
ax.LineWidth = 2;
grid on
ax.FontWeight = 'bold';
set(ax, 'LineWidth', 3);
ax = gca;
ax.FontSize =20;  
set(gca, 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











 