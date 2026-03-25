function [K_opt,z_opt,elapsed_time,pars,obj,varargout]=tpddae_stabopt_static(rddae,M,type,varargin)

%%
% This function returns the optimal time-periodic feedback gain, the objective
% function, and the elapsed time used to solve the optimization problem via HANSO.
%
% Required Inputs: 
% rddae: A struct describing a time-periodic delay differential equation of the retarded type.
%
% M: The degree of the polynomial approximating the solution in each subinterval. 
%
% initial: Defines how the initial starting point of HANSO should be chosen.
%
% type: The type of the designed controller, which can be:
%   - 'TC': Time-constant
%   - 'QR': Quadratic form
%   - 'FR': Fourier form
%
% weight: The weight function in the optimization problem for the 'QR' method. 
% 
% Depending on the type, you must assign the required optimization
% parameters: 
%   - If the feedback gain type is 'QR', you must assign a value for 'lambda_1'.
%   - If the feedback type is 'FR' (Fourier), you must assign the number of Fourier terms via 'n_Fr',
%     which must be an odd number.
%
% Outputs:
% K_opt: The optimal feedback gain.
% z_opt: The spectral radius of the monodromy operator using K_opt.
%
% Optional Inputs: 
% method: The method for solving the general eigenvalue problem, which is
% either 'eigs' or 'eig'. 
%
% Delta: The shortest length of the subinterval for approximating
% the solution. If there is a breakpoint in the interval [0,T],
% it has priority in determining the subintervals. If not specified, 
% the length of the subinterval will be equal to T. 
%
% normtol: Termination tolerance for the norm of the smallest vector 
% in the convex hull of a set of saved gradients.
%
% n_start: The maximum number of iterations for each BFGS starting point. 
% This must be greater than 10.
%
% maxit_gradsam: The maximum number of iterations for gradient sampling.
%
% initial: Columns represent one or more starting vectors of variables used to 
% initialize the BFGS phase of the algorithm.

%%
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
pars=pars_FM(rddae,M,type,opts);
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
[K,obj]=hanso(pars,options);
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
[z_opt,~]=FM_computation(rddae,M,'gain',gain,'method',FLAG,'slack',pars.slack);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if (abs(z_opt)>1)
         warning('The obtained optimal feedback gain is not stabilizing')
     end

    elapsed_time=toc;
    
    
    Delta=pars.Delta;
    K_opt=controller_fun(rddae,M,gain,Delta); 
end




