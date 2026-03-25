function [abs_obj,g_obj]=opt_func_UPO(v_gain,pars)
%%% v_gain is a the feedback matrix in a vector of size r*p*(M)*NT. Thus, we need to reshape it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rddae=pars.rddae;
M=pars.M;
type=pars.type;
FLAG=pars.FLAG;
B=rddae.B;
r=size(B,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain_=reshape(v_gain,r,[]);
gain.type=type;
gain.K=gain_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if type=='TC'
    
     [z,g_z]=FM_computation_UPO(rddae,M,'gain',gain,'method',FLAG,'slack',pars.slack);
    abs_obj=abs(z); %%% Calculating the value of the objective function
    g_z=g_z(:);
    g_obj=g_z;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
elseif type=='FR'

    gain.n_FR=pars.n_FR;
    [z,g_z]=FM_computation_UPO(rddae,M,'gain',gain,'method',FLAG,'slack',pars.slack);
    abs_obj=abs(z);%%% Calculating the value of the objective function
    g_z=g_z(:);
    g_obj=g_z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif type=='QR'


    lambda1=pars.lambda1; 
    [z,g_z]=FM_computation_UPO(rddae,M,'gain',gain,'method',FLAG,'slack',pars.slack);
    LAMBDA_1=pars.LAMBDA_1;
    gain_=gain_(:);
    gain_2=gain_;
    g_reg1=lambda1*(LAMBDA_1+LAMBDA_1')*gain_2;
    obj=trace(gain_2'*LAMBDA_1*gain_2);
    abs_obj=abs(z)+lambda1*obj;
    g_reg=g_reg1(:);
    g_z=g_z(:);
    g_obj=g_z+g_reg;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
