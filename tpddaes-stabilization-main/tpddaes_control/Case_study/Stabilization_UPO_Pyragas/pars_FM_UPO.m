function pars=pars_FM_UPO(rddae,M,type,options)
% This file returns strcture pars which contains the information required 
% for running the optimization problem with HANSO. 
% Input Paramters:
% rddae: A struct containing the Time-periodic Delay Differential Algebraic Equation of retarded type
% M: degree of the polynomial approximating the solution 
% type: design type of the time-periodic feedback gain. It can be either
% 'TC' standing for time-constant, 'QR' standing for Quadaratic form, and
% finally, 'FR' standing for Fourier form.
% options: an optional struct which contains:
% options.method: method of solving the general eigenvalue problem which
% can be either 'eigs' or 'eig'.
% options.Delta: the shortest length of the subinterval for apporixmating
% the solution. Note that if there is a break point in the interval [0,T],
% it has prioirty in choosing the subintervals. If not indicated the length
% of the subinterval will be equal to T. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B=rddae.B;
C=rddae.C;
r=size(B,2);
p=size(C,1);
pars.rddae=rddae;
pars.type=type;
pars.M=M;
T=rddae.T;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(options,'method')
FLAG=options.method;
else 
    FLAG='eig';
end
pars.FLAG=FLAG;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(options,'Delta')
  Delta=options.Delta;
else 
    Delta=T;
end
pars.Delta=Delta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if type =='FR'
        
          n_FR=options.n_FR;
 else 
   n_FR=[];
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slack_out=slack_function(Delta,M,rddae,type,n_FR);
pars.slack=slack_out;
NT=slack_out.NT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pars.fgname = 'opt_func_UPO';
if type=='TC'
    pars.nvar=r*p;
elseif type=='FR'
    n_FR=options.n_FR;
    pars.nvar=r*p*n_FR;
    pars.n_FR=n_FR;
elseif type=='QR'
    M_=M+1;
    pars.nvar=r*p*(M_)*NT-r*p*(NT);
    pars.lambda1=options.lambda1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(options,'weight')
        weight=options.weight;
    else 
        weight=@(t) 1;
    end
    [LAMBDA_1,slack]=QR_func(rddae,M,weight,slack_out);
    pars.LAMBDA_1=LAMBDA_1;
    pars.slack_mat=slack;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end