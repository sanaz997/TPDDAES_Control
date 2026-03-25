function slack=slack_function(Delta,M,rddae,type,n_FR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=rddae.A;
tau_A=rddae.tau_A;
points_d=rddae.points_d;
tau_s={@(t) 0};
tau=[tau_A ;tau_s];
m=length(tau);
d=size(A{1}(0),1); 
T=rddae.T;
M_=M+1;
E=rddae.E;
C=rddae.C;
p=size(C,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intr_0=[0 points_d T] ;
clear interval_d
s=1;
for i=1:length(intr_0)-1  
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Largest time delay:
tau_star=zeros(1,m);
for i=1:m
    f=@(t) t-tau{i}(t);
    opts.TolX = 1.0E-13;
    [~, tau_star(i)] = fminbnd(f, 0, T, opts);
end
tau_max=-min(tau_star);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Choosing the Chebyshev points in the interval [0,T]
alphaR=(-cos((1:M)*pi/(M)));

intR=[pointsR(1:end-1)' pointsR(2:end)'];
thetaR=cell(1,NT);
for i=1:size(intR,1)
    thetaR{i}=0.5*(intR(i,1)+intR(i,2))+alphaR/2*(intR(i,2)-intR(i,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Choosing the Chebyshev points in the interval [-\tau_max,0]
alphaL=(-cos((0:M)*pi/(M)));
cc=floor(tau_max/T);       % largest integer s.t. cc*T<=tau_max
Kl=cc*NT;                   % number of subintervals of [-tau_max,0] IF cc*T=tau_max
dist=tau_max-cc*T;
negpoints=points-T;
if dist>0                   % if tau_max>cc*T, then we look for the last discontinuity point P such that -tau_max<-cc*T-P
    k=1;
    while (k<=length(points)) && dist>-negpoints(end-k+1) 
        k=k+1;
    end
    Kl=Kl+k;
end


if cc*T~=tau_max
    pointsL=[-tau_max -cc*T+negpoints(end-k+2:end) -cc*T];
else
    pointsL=-tau_max;
end

if cc>0
    l=cc;
    while l>0
    pointsL=[pointsL -(l-1)*T+negpoints -(l-1)*T];
    l=l-1;
    end
end
intL=[pointsL(1:end-1)' pointsL(2:end)'];

thetaL=cell(1,Kl);
for i=1:size(intL,1)
    thetaL{i}=0.5*(intL(i,1)+intL(i,2))+alphaL/2*(intL(i,2)-intL(i,1)); % every vector in thetaL defines a grid of M Chebyshev-distributed points on each subinterval of [-tau_m,0]
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:m-1

    for l=1:NT
        for j=1:M_-1                                      
                 Dvalues{i,j,l}=A{i}(thetaR{l}(j));
        
        end
      
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if type=='FR'
%     n_FR=gain.n_FR;
thetaR_matrix = cell2mat(thetaR); % Convert cell array to matrix
thetaR_tot = reshape(thetaR_matrix, [], 1);

FR_matrix=zeros(n_FR,length(thetaR_tot));

FR_matrix(1,:)=1;
     for s=2:2:n_FR
        FR_matrix(s,:)=cos(2*pi*s/2*thetaR_tot(:)/T);
        FR_matrix(s+1,:)=sin(2*pi*s/2*thetaR_tot(:)/T);
     end
    
else 
FR_matrix=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here we define collocation constraints for the definition of eigenfunction: p(t+T)=mu*p(t) for t in [-tau_m,0];
Bp=zeros(M_*(Kl+NT));
for i=1:M_
    for k=1:Kl
        Bp((k-1)*M_+1:k*M_,M_*(k-1)+i)=chebT(alphaL,i-1)';
    end
end
Bp=kron(Bp,eye(d));         % this defines the left-hand side of the equation

Ap=zeros(M_*(Kl+NT));

int=[intL;intR];
for i=1:M_           
    for k=NT:Kl+NT-1
        Ap(M_*(k-NT)+1:M_*(k-NT+1),M_*k+i)=chebT(2/(int(k+1,2)-int(k+1,1))*(thetaL{k-NT+1}+T-0.5*(int(k+1,2)+int(k+1,1))),i-1)';
    end
end
Ap=kron(Ap,eye(d));         % this defines the right-hand side of the equation
%% Here we impose continuity in 0 and in any breaking point of the subinterval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clear rowL rowR
for i=1:M+1
    rowL(i)=chebT(1,i-1);
end
for i=1:M+1
    rowR(i)=-chebT(-1,i-1);
end
rowL=kron(rowL,eye(d));
rowR=kron(rowR,eye(d));

Bp(M_*d*Kl+1:M_*d*Kl+d,M_*d*(Kl-1)+1:M_*d*(Kl-1)+M_*d)=rowL;        %rightmost subinterval in [-tau,0] 
Bp(M_*d*Kl+1:M_*d*Kl+d,M_*Kl*d+1:M_*Kl*d+M_*d)=rowR;      %leftmost subinterval in [0,T]

for i=1:NT-1
Bp(M_*d*Kl+d+(i-1)*d+1:M_*d*Kl+d+i*d,M_*d*Kl+d*(i-1)*M_+1:M_*d*Kl+d*i*M_)=rowL;
Bp(M_*d*Kl+d+(i-1)*d+1:M_*d*Kl+d+i*d,M_*d*Kl+d*i*M_+1:M_*d*Kl+d*(i+1)*M_)=rowR;
end

%% Here we define collocation constraints for the TPDDE
P0=zeros((M_-1)*NT,M_*Kl+M_*NT);
for k=1:NT
for i=2:M_

P0((M_-1)*(k-1)+1:(M_-1)*k,Kl*M_+(k-1)*M_+i)=2/(intR(k,2)-intR(k,1))*(i-1)*chebU(alphaR,i-2)';  
end
end

P0=kron(P0,E);     % In P0 we store all coefficient matrices multiplying the derivative p'(t) (this is why we include Chebyshev polynomials of the second type)


for s=1:NT
for k=1:M_-1
    for l=1:m
        c=thetaR{s}(k)-tau{l}(thetaR{s}(k));
        for j=1:NT+Kl
        if c>int(j,1) && c<=int(j,2)
            
            K_index(s,k,l)=j;
            break
        end
        end
    end
end
end




%% new part of the code
P=cell(1,m-1);


for k=1:m-1
    P{k}=zeros(d*(M_-1)*NT,d*M_*Kl+d*M_*NT);
  
    for l=1:NT
    for j=1:M_-1
        cc=K_index(l,j,k);
        p1=int(cc,1);
        p2=int(cc,2);
         for i=1:M_


% Dvalues{k,j,l}

        P{k}((l-1)*d*(M_-1)+(j-1)*d+1:(l-1)*d*(M_-1)+j*d,d*M_*(cc-1)+(i-1)*d+1:d*M_*(cc-1)+d*i)=Dvalues{k,j,l}*chebT(2/(p2-p1)*(thetaR{l}(j)-tau{k}(thetaR{l}(j))-0.5*(p1+p2)),i-1);
      
         end
    end
    end
end                % In P{k} we store all the coefficient matrices multiplying p(t-tau{k}(t))  





finalP=P0;

for i=1:m-1
    finalP=finalP-P{i};
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 P2=zeros(p*(M)*NT,d*M_*Kl+d*M_*NT);
    for l=1:NT
    for j=1:M_-1
        cc=l+Kl;
        p1=int(cc,1);
        p2=int(cc,2);
        for i=1:M_
 P2((l-1)*p*(M_-1)+(j-1)*p+1:(l-1)*p*(M_-1)+j*p,d*M_*(cc-1)+ ...
            (i-1)*d+1:d*M_*(cc-1)+d*i)=C*chebT(2/(p2-p1)*(thetaR{l}(j)-0.5*(p1+p2)),i-1);       
        end
    end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_K2=zeros(M_*p*(Kl+NT),(M_)*d*(Kl+NT));
M_K2(M_*Kl*p+NT*p+1:end,:)=-P2;
Matrix=M_K2((M_)*Kl*p+NT*p+1:end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% making the output
slack= struct('Ap', Ap, 'Bp', Bp, 'NT', NT,'Kl',Kl,'thetaR',thetaR,'K_index',K_index,'int',int,'finalP',finalP,'Matrix',Matrix,'FR_matrix',FR_matrix,'Delta',Delta,'intR',intR);
