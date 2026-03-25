function plot_K(K_fun)
% This function plots the feedback gain value with respect to time over the interval [0,T]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if K_fun.type=='QR' %%% Feedback gain designed based on total variation method

points_d=K_fun.points_d;
K_val=K_fun.K;
[r,p]=size(K_val.K_1);
Delta=K_fun.Delta;
NT=K_fun.NT;
M=K_fun.M;
Theta=K_fun.theta;
fieldName2 = sprintf('theta_%d', M*NT);
T=K_fun.theta.(fieldName2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% We have to make slack_mat:
intr_0=[0 points_d T] ;
clear interval_d
s=1;
for i=1:length(intr_0)-1
if int32(intr_0(i+1)-intr_0(i))<=Delta
    intr_d(1,s)=intr_0(i);
   s=s+1;
else 
    domain=(intr_0(i+1)-intr_0(i));
    num=floor(domain/Delta);
   intr_d(1,s:(s+num))=intr_0(i)+(0:1:num)*Delta;
s=s+(num+1);
end
end

if (intr_d(end)~=T)
intr_d(end+1)=T;
end
points= intr_d(1,2:end-1);
Dis2=[0 points T];
pointsR=Dis2;
NT=length(Dis2)-1;
alphaR=(-cos((1:M)*pi/(M)));
intR=[pointsR(1:end-1)' pointsR(2:end)'];
thetaR=cell(1,NT);
for i=1:size(intR,1)
    thetaR{i}=0.5*(intR(i,1)+intR(i,2))+alphaR/2*(intR(i,2)-intR(i,1));
end








M_=M+1;



for s=1:NT
for k=1:M 
        c=thetaR{s}(k);
        for j=1:NT
        if c>intR(j,1) && c<=intR(j,2)
            K_index(s,k,1)=j;
            break
        end
    end
end
end       
%%% Q: represents the polynomial approximating the behvaior of the feedbackgain:
%%%% writing the action of interpolation operator
rp=1;
for l=1:NT
    for j=1:M_-1
        cc=K_index(l,j,1);
        p1=intR(cc,1);
        p2=intR(cc,2);
        for i=1:M_

      Q{k}((l-1)*rp*(M_-1)+(j-1)*rp+1:(l-1)*rp*(M_-1)+ ...
          j*rp,rp*M_*(cc-1)+(i-1)*rp+1:rp*M_*(cc-1)+rp*i)=chebT(2/(p2-p1)*(thetaR{l}(j)-0.5*(p1+p2)),i-1);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
finalQ=0;
for i=1:length(Q)
finalQ=Q{i};
end
clear slack 
slack=zeros((M+1)*(NT));
slack(NT+1:end,:)=(finalQ);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:NT-1
slack((i-1)+1:i,(i-1)*M_+1:i*M_)=chebT(1,0:M); %rightmost point of subinterval number i
slack((i-1)*rp+1:i*rp,rp*i*M_+1:rp*(i+1)*M_)=-chebT(-1,0:M); %%leftmost point of subinterval number i
end

slack(NT,1:M_)=-chebT(-1,0:M) ;
slack(NT,(NT-1)*M_+1:NT*M_)=slack(NT,(NT-1)*M_+1:NT*M_)+chebT(1,0:M) ; %%%% T
slack_mat=kron(slack,eye(r*p));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure();

i=1;
time=0;

for j=0:T/1000:T
   time(i)=j;
    gain2=K_continous(K_fun,j,intR,slack_mat);
 gain(1:r*p,i)=gain2(:); %%% stroing value of the feedback gain over the time period of [0,T]
 i=i+1;
end

for k=1:r*p
 plot(time,gain(k,:),'LineWidth',5)
     hold on;
end

grid on;
xlabel('Time-interval','FontSize',16, 'FontWeight','normal')
%% title('Design based on the polynomial description of the feedback gain over a time-period','FontSize',16)

title('QR')

%%% legend:
s=1;
for j = 1:p
    for i=1:r
   legendInfo{s} =['\boldmath $\mathcal{K}_{' num2str(i) , num2str(j) '}$'];
   s=s+1;
    end
end
lgd=legend(legendInfo, 'interpreter', 'latex');

fontsize(lgd,25,'points');

xlim([0 T])
ax = gca;
ax.FontSize = 25;  % Font Size of 15
xlabel('t(s)','FontSize',25,'FontWeight', 'bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif K_fun.type=='TC' %%% Feedback gain designed based on truncated Fourier method
T=K_fun.T;

[r,p]=size(K_fun.K);
gain=K_fun.K;



x=1;
for j=0:T/100:T    

val(1:r*p,x)=gain(:); %%% Value of the feedback gain
vec_x(x,1)=j; %%% Time vector
x=x+1;
end

% 
 for k=1:r*p
     plot(vec_x,val(k,:),'LineWidth',5)
     hold on
 end


grid on;
xlabel('Time-interval','FontSize',16, 'FontWeight','bold')
ylabel('\boldmath  $\mathcal{K}(t)$ ','FontSize',16, 'interpreter', 'latex')
title('TI' ,'FontSize',16,'FontWeight','bold')


%%% legend:
s=1;
for j = 1:p
    for i=1:r
    legendInfo{s} =['\boldmath $\mathcal{K}_{' num2str(i) , num2str(j) '}$'];
     s=s+1;
    end
end
lgd=legend(legendInfo, 'interpreter', 'latex');
fontsize(lgd,25,'points');
ax = gca;
ax.FontSize = 25;  

xlim([0 T])
xlabel('t(s)','FontSize',25,'FontWeight', 'bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif K_fun.type=='FR' %%% Feedback gain designed based on truncated Fourier method
T=K_fun.T;
n_FR=K_fun.n_FR;
f_gain = cell(n_FR,1);
nn=(n_FR-1)/2;
f_gain{1}=@(t) 1;

[r,p]=size(K_fun.K_1);

if nn~=0
for s=2:2:n_FR

        f_gain{s,1}=@(t) cos((2*pi*s/2*t)/T);
        f_gain{s+1,1}=@(t) sin((2*pi*s/2*t)/T);
end
end




% figure();
for i = 1:n_FR
     fieldName = sprintf('K_%d', i);
    gain_(1:r,p*(i-1)+1:p*i)=K_fun.(fieldName);  %%% storing feedback gain from the structure array
end





x=1;
for j=0:T/100:T    
sum_val=0;
for s=1:n_FR

        eval_Fr=(gain_(1:r,p*(s-1)+1:p*s))*f_gain{s}(j);
        sum_val=sum_val+eval_Fr; %%% stroing value of the feedback gain over the time period of [0,T]
end

val(1:r*p,x)=sum_val(:); %%% Value of the feedback gain
vec_x(x,1)=j; %%% Time vector
x=x+1;
end

% 
 for k=1:r*p
     plot(vec_x,val(k,:),'LineWidth',5)
     hold on
 end

grid on;
xlabel('Time-interval','FontSize',16, 'FontWeight','bold')
ylabel(' \boldmath $\mathcal{K}(t)$ ','FontSize',16, 'interpreter', 'latex')

title('FR')

%%% legend:
s=1;
for j = 1:p
    for i=1:r
legendInfo{s} =['\boldmath $\mathcal{K}_{' num2str(i) , num2str(j) '}$'];
    s=s+1;
    end
end

lgd=legend(legendInfo, 'interpreter', 'latex');
fontsize(lgd,25,'points');
ax = gca;
ax.FontSize = 25;  
xlabel('t(s)','FontSize',25,'FontWeight', 'bold')


end



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
xlim([0 T])


















 