function K_fun=controller_fun(rdde,M,gain,Delta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if gain.type=='FR'
type=gain.type;
K_fun.type=type;

B=rdde.B;
C=rdde.C;
[r]=size(B,2);
[p]=size(C,1);
n_FR=gain.n_FR;
gain2=gain.K;
T=rdde.T;
K_fun.T=T;
if type=='FR'
    K_fun.n_FR=n_FR;

    for i = 1:n_FR
     fieldName = sprintf('K_%d', i);
    K_fun.(fieldName) =(gain2(1:r,p*(i-1)+1:p*i));
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_fun.type=='FR' %%% Feedback gain designed based on truncated Fourier method
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








elseif gain.type=='QR'
K_TV=gain.K;

B=rdde.B;
C=rdde.C;
[r]=size(B,2);
[p]=size(C,1);


T=rdde.T;
points_d=rdde.points_d;

intr_0=[0 points_d T] ;
clear interval_d
s=1;
for i=1:length(intr_0)-1  
   %(intr_0(i+1)-intr_0(i))
if (intr_0(i+1)-intr_0(i) )==Delta || (intr_0(i+1)-intr_0(i) )<Delta
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


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 alphaR=(-cos((1:M)*pi/(M)));
intR=[pointsR(1:end-1)' pointsR(2:end)'];
thetaR=cell(1,NT);

for i=1:size(intR,1)
    thetaR{i}=0.5*(intR(i,1)+intR(i,2))+alphaR/2*(intR(i,2)-intR(i,1));
end
thetaR{end}(end)=T; 

M_=M+1;
for i=1:NT
theta_tot(1+(i-1)*(M_-1):i*(M_-1))=thetaR{i}(1:M_-1);
end


gain2=reshape(K_TV,r,[]);
for i = 1:length(theta_tot)
     fieldName = sprintf('K_%d', i);
   fieldName2 = sprintf('theta_%d', i);
    K_fun.(fieldName) =(gain2(1:r,p*(i-1)+1:p*i));
    K_fun.(fieldName2)=theta_tot(i);
end

K_fun.NT=NT;
K_fun.M=M;
K_fun.points_d=rdde.points_d;
K_fun.Delta=Delta;
K_fun.type='QR';
K_fun.len=length(theta_tot);

%%% 
K=struct();

for i = 1:length(theta_tot)
     fieldName = sprintf('K_%d', i);
   fieldName2 = sprintf('theta_%d', i);
   K.(fieldName) =(gain2(1:r,p*(i-1)+1:p*i));
   theta.(fieldName2)=theta_tot(i);
end

K_fun=struct();
K_fun.K=K;
K_fun.theta=theta;

K_fun.NT=NT;
K_fun.M=M;
K_fun.points_d=rdde.points_d;
K_fun.Delta=Delta;
K_fun.type='QR';
K_fun.len=length(theta_tot);
K_fun.T=T;





elseif gain.type=='TC'

type=gain.type;
K_fun.type=type;

B=rdde.B;
C=rdde.C;
[r]=size(B,2);
[p]=size(C,1);
gain2=gain.K;
T=rdde.T;
K_fun.T=T;
gain2=reshape(gain2,r,[]);
K_fun.K =(gain2(1:r,1:p));



end


