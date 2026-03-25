function Cheb_vec=func_coef(slack_mat,K_fun)
NT=K_fun.NT; %%% number of subintervals in [0,T]
M=K_fun.M;
Dyn_flag=isfield(K_fun,'Ac');
Delay_flag=isfield(K_fun,'tau_1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Dyn_flag==0 && Delay_flag==0
%%% extracting the values of the feedback gain from the structure array
    for i = 1:(M)*NT
         
    
        fieldName = sprintf('K_%d', i);
        val=K_fun.K;
        [r,p]=size(val.K_1); %%% size of matrix K
        
        KK=val.(fieldName);
        K_TV(r*p*(i-1)+1:r*p*i,1)=KK(:);
    
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif Dyn_flag==1
    
    theta_vec=K_fun.theta;
    for i=1:length(fieldnames(theta_vec))
    fieldName = sprintf('K_%d', i);
    Ac_mat=K_fun.Ac;
    Bc_mat=K_fun.Bc;
    Cc_mat=K_fun.Cc;
    Dc_mat=K_fun.Dc;
    a=Ac_mat.(fieldName);
    b=Bc_mat.(fieldName);
    c=Cc_mat.(fieldName);
    d=Dc_mat.(fieldName);
    [q,p]=size(b);
    [r,q]=size(c);
    r=1;
    fieldName2 = sprintf('theta_%d', i);  
    thetaRR(i)= theta_vec.(fieldName2);
    col=q+p;
    row=q+r;
    K_TV(1:row,col*(i-1)+1:col*i)=[a b;c d];
    end

elseif Delay_flag==1
    K_TV=[];


 theta_vec=K_fun.theta;

 for i=1:length(fieldnames(theta_vec))

 for j=1:K_fun.len_tau
        fieldName3= sprintf('tau_%d', j); 
        
        slack=K_fun.(fieldName3);
        
        
        fieldName = sprintf('K_%d', i);
        a=slack.(fieldName);
        [r,p]=size(a);
        r_=r;p_=p;
        K_TV=[K_TV(:);a(:)];
        
        fieldName2 = sprintf('theta_%d', i);  
        thetaRR(i)= theta_vec.(fieldName2);
        col=p*K_fun.len_tau;
        row=r;
end

end
end









if Dyn_flag==1

r=row;
p=col;
end



if Delay_flag==1

r=row;
p=col;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B0=zeros(r*p*NT,1);
B1=K_TV;
B_=[B0;B1];
Cheb_vec=inv(slack_mat)*B_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end