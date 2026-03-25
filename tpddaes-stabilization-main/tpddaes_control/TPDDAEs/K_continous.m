function [output,NT]=K_continous(K_fun,time,intR,slack_mat)
 NT=K_fun.NT; %%% number of subintervals in [0,T]
M=K_fun.M;
K_val=K_fun.K;
[r,p]=size(K_val.K_1);
cheb_vec=func_coef(slack_mat,K_fun);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% finding variable time falls in which subinterval and
%%%% storing the upper and lower bounds in p2 and p1, respectively
for s=1:NT
        c=time;
        for j=1:NT
        if c>=intR(j,1) && c<=intR(j,2)
p2=intR(j,2);
    p1=intR(j,1);
index=j;
            break
        end
        end
end       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cheb_poly=zeros(1,(M+1)*NT);
cheb_poly(1,(index-1)*(M+1)+1:(M+1)*index)=chebT((2/(p2-p1))*(time-0.5*(p1+p2)),[0:1:M]);
cheb_poly=kron(cheb_poly,eye(r*p));
v_gain=cheb_poly*cheb_vec;
output=reshape(v_gain,r,[]);
