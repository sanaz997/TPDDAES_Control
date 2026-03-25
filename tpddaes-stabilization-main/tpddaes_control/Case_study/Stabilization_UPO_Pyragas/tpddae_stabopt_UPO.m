function [K_opt,z_opt,elapsed_time,pars,obj,initial,varargout]=tpddae_stabopt_UPO(rddae,M,type,varargin)

% This function returns the optimal time-periodic feedback gain, objective
% function and the elapsed time used to solve the optimization problem via HANSO.
%
% Required Inputs: 
% rddae: a struct describing a time-periodic delay differential equation of retarded type.
%
% M: is the degree of polyunomial approximating the solution in each subinterval 
% initial defines how to initial starting point of HANSO should be chosen.
%
% type: type of the designed controller which can be 'TC', standing for
% time-constnat, 'QR' standing for Quadratic, or 'FR' standing for Fourier
% form.
%
% weight: stands for the weight function in the optimization problem for 'QR' method. 
% 
% Depending on type, you have to assign the required optimization
% paramters: 
% If the type of the feedback gain is 'QR', you have to assign value of
% 'lambda_1'. If the type of the feedback is Fourier you have to assign the
% number of Fourier via 'n_Fr' terms which has to be odd.
%
% Outputs:
% K_opt: optimal feedback gain
% z_opt: spectral radius of the monodromy operator using K_opt
%
%
% You can also assign optional inputs: 
% method: the method of solving the general eigenvalue problem, which is
% either 'eigs' or 'eig'. 
%
% Delta: the shortest length of the subinterval for apporixmating
% the solution. Note that if there is a break point in the interval [0,T],
% it has prioirty in choosing the subintervals. If not indicated the length
% of the subinterval will be equal to T. 
%
% normtol: termination tolerance for norm of smallest vector in convex hull of a set of saved gradients:
% n_start: max number of iterations for each BFGS starting point. You can
% only assign more than 10.
%
% maxit_gradsam: max number of iterations for gradient sampling.
%
% initial:columns are one or more starting vector of variables used to 
% intialize the BFGS phase of the algorithm
%


tic
Names = { 'initial', 'Delta','method','lambda1','n_FR','nstart','maxit_gradsam' ,'maxit','normtol','weight'};
m= length(Names);
opts = [];
i = 1;
num=nargin-3;
while i <= num
  arg = varargin{i};
    if ~isempty(arg)          
        flag=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:m
            if any(strcmp(arg,Names{j})) 
                val=varargin{i+1};   
                flag=0;
            else
                val=[];
            end
            if ~isempty(val)
            opts.(Names{j}) =  val;
            end
        end
            if flag==1
            error(['No such field available as: ', num2str(arg)]);
            end

    end
    i = i + 2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('options')
opts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pars=pars_FM_UPO(rddae,M,type,opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(opts,'maxit_gradsamp')
    options.maxit_gradsamp =opts.maxit_gradsamp;
else 
    options.maxit_gradsamp = 0;
end
if isfield(opts,'initial')
    options.x0 =[opts.initial randn(pars.nvar, opts.nstart-1)];
end
if isfield(opts,'normtol')
    options.normtol= opts.normtol;
end

if isfield(opts,'nstart')
    options.nstart=opts.nstart;
end
if isfield(opts,'maxit')
    options.maxit=opts.maxit;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[K,obj,initial]=hanso(pars,options);
if type=='QR'
    LAMBDA_1=pars.LAMBDA_1;
    reg_term=trace(K'*LAMBDA_1*K);
    varargout{1}=reg_term;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=rddae.B;
r=size(B,2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if type=='FR'
    gain.n_FR=pars.n_FR;
end
gain_=reshape(K,r,[]);
gain.K=gain_;
gain.type=type;
FLAG=pars.FLAG;
[z_opt,~]=FM_computation_UPO(rddae,M,'gain',gain,'method',FLAG,'slack',pars.slack);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if (abs(z_opt)>1)
         warning('The obtained optimal feedback gain is not stabilizing')
     end

elapsed_time=toc;
Delta=pars.Delta;
K_opt=controller_fun(rddae,M,gain,Delta); 

end




