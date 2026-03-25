function  [rddae_2,gain_final,z,Time,pars,gain,reg]=tpddae_stabopt_dynamic(rddae,q,M,type,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=rddae.A;
B=rddae.B;
E=rddae.E;
C=rddae.C;
tau_A=rddae.tau_A;
T=rddae.T;
points_d=rddae.points_d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is for desisning a dynamic feedback controller:
r=size(B,2);
[p,d]=size(C);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E_=[E zeros(d,r) zeros(d,q);
   zeros(q,d) zeros(q,r) eye(q,q) ;
   zeros(r,d) zeros(r,r) zeros(r,q)];


l_A=length(A);
A_=cell(l_A,1);
tau_A_=cell(l_A,1);

A_{1}=@(t) [A{1}(t) B zeros(d,q);
            zeros(q,d) zeros(q,r) zeros(q,q) ;
            zeros(r,d) -eye(r,r) zeros(r,q)  ];

tau_A_{1}=@(t) 0;


for i=2:l_A
    
    A_{i,1}=@(t) [A{i}(t) zeros(d,r) zeros(d,q);
                 zeros(q,d) zeros(q,r) zeros(q,q) ;
                 zeros(r,d) zeros(r,r) zeros(r,q) ];


    tau_A_{i,1}=@(t) tau_A{i}(t);
     
end


C_= [zeros(q,d)  zeros(q,r) eye(q,q);
     C   zeros(p,r) zeros(p,q)];

B_= [zeros(d,q) zeros(d,r);
    eye(q,q) zeros(q,r);
    zeros(r,q) eye(r,r)];


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% q: size feedback gain;
R=null(E_');
Q=null(E_);

 if(norm(R'*B_)~=0 && norm(C_*Q)~=0)
     warning('Assumption 2 is not satisfied.')
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rddae_2=tpddae_create(A_,B_,C_,E_,tau_A_,points_d,T);
[gain,z,Time,pars,reg]=tpddae_stabopt_static(rddae_2,M,type,varargin{:})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AC=struct();
BC=struct();
CC=struct();
DC=struct();
gain_final=struct();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if gain.type=='FR'

    for i=1:gain.n_FR
    
        fieldName = sprintf('K_%d', i);
        val=gain.(fieldName);
        AC.(fieldName) =val(1:q,1:q);
        BC.(fieldName) =val(1:q,q+1:end);
        CC.(fieldName) =val(q+1:end,1:q);
        DC.(fieldName) =val(q+1:end,q+1:end);
    end

gain_final.n_FR=gain.n_FR;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif gain.type=='TC'

    val=gain.K;
    AC =val(1:q,1:q);
    BC =val(1:q,q+1:end);
    CC =val(q+1:end,1:q);
    DC =val(q+1:end,q+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif gain.type=='QR'
    theta=struct();

    for i=1:M*gain.NT
    
        fieldName = sprintf('K_%d', i);
        slack=gain.K;
        val=slack.(fieldName);
        AC.(fieldName) =val(1:q,1:q);
        BC.(fieldName) =val(1:q,q+1:end);
        CC.(fieldName) =val(q+1:end,1:q);
        DC.(fieldName) =val(q+1:end,q+1:end);
        fieldName2 = sprintf('theta_%d', i);  
        theta_vec=gain.theta;
        theta.(fieldName2)=theta_vec.(fieldName2)';

    end


gain_final.M=gain.M;
gain_final.NT=gain.NT;
gain_final.theta=theta;

 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gain_final.type=gain.type;
gain_final.Ac=AC;
gain_final.Bc=BC;
gain_final.Cc=CC;
gain_final.Dc=DC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%% This function plots the optimal feedback gain with respect to t\in [0,T].
if type=='QR'
    plot_K_dynamic(rddae_2,q,gain,pars.slack, pars.slack_mat) 
else
    plot_K_dynamic(rddae_2,q,gain,pars.slack)
end


end
















