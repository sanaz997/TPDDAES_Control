function pars=pars_FM(rddae,M,type,options)

% This file returns the structure 'pars', which contains all the necessary 
% information required to run the optimization problem using HANSO.
%
% Input Parameters:
% -----------------
% rddae:    A struct representing the Time-Periodic Delay Differential-Algebraic 
%           Equations (TPDDAEs) of the retarded type.
%
% M:        The degree of the polynomial used to approximate the solution.
%
% type:     The design type of the time-periodic feedback gain. It can be:
%              - 'TC' : Time-Constant
%              - 'QR' : Quadratic Form
%              - 'FR' : Fourier Form
%
% options:  (Optional) A struct containing additional configuration parameters:
%              - options.method: The method for solving the general eigenvalue 
%                                problem, either 'eigs' or 'eig'.
%
%              - options.Delta: The shortest length of the subinterval used to 
%                                approximate the solution. If there is a breakpoint 
%                                within the interval [0, T], it takes priority 
%                                when selecting subintervals. If not specified, 
%                                the subinterval length will default to T.

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
pars.fgname = 'opt_func';
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