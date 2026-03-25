function [LAMBDA_1,slack]=QR_func(rddae,M,weight,slack_out)
M_=M+1;
B=rddae.B;
C=rddae.C;
r=size(B,2); 
p=size(C,1); 
clear P;clear B_;clear v_ ;clear v;clear K_TV finalQ Q alphaR slack intr_d F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[NT,intR]=deal(slack_out(1).NT,slack_out(1).intR);
thetaR = cell(1, NT);
[thetaR{:}]=deal(slack_out(1:NT).thetaR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% finding the correct interval for thetaR{s}(k) and storing it in K_index(s,k,1)
for s=1:NT
for k=1:M 
        c=thetaR{s}(k);
        for j=1:NT
        if c>intR(j,1) && c<=intR(j,2)
            K_index(s,k,1)=j;
            break
        end
    end
end
end       
%%% Q: represents the polynomial approximating the behvaior of the feedback gain:
%%%% writing the action of interpolation operator
rp=1;
for l=1:NT
    for j=1:M_-1
        cc=K_index(l,j,1);
        p1=intR(cc,1);
        p2=intR(cc,2);
        for i=1:M_

      Q{k}((l-1)*rp*(M_-1)+(j-1)*rp+1:(l-1)*rp*(M_-1)+ ...
          j*rp,rp*M_*(cc-1)+(i-1)*rp+1:rp*M_*(cc-1)+rp*i)=chebT(2/(p2-p1)*(thetaR{l}(j)-0.5*(p1+p2)),i-1);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
finalQ=0;
for i=1:length(Q)
finalQ=Q{i};
end
clear slack 
slack=zeros((M+1)*(NT));
slack(NT+1:end,:)=(finalQ);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:NT-1
slack((i-1)+1:i,(i-1)*M_+1:i*M_)=chebT(1,0:M); %rightmost point of subinterval number i
slack((i-1)*rp+1:i*rp,rp*i*M_+1:rp*(i+1)*M_)=-chebT(-1,0:M); %%leftmost point of subinterval number i
end

slack(NT,1:M_)=-chebT(-1,0:M) ;
slack(NT,(NT-1)*M_+1:NT*M_)=slack(NT,(NT-1)*M_+1:NT*M_)+chebT(1,0:M) ; %%%% T
slack=kron(slack,eye(r*p));
inv_slack=inv(slack);
inv_slack(:,1:NT*r*p)=[];Matrix_N=inv_slack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%% Objective function limiting the total variation of the feedback gain:
 F3= integral(@(t) d_cheb2(t,M,weight),-1,+1,'ArrayValued', true);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%M=
for i=1:NT

UTU((i-1)*M_+1:i*M_,(i-1)*M_+1:i*M_)=(intR(i,2)-intR(i,1))/2*F3; %%% size MN_T\times MN_T
end
UTU=kron(UTU,eye(r*p));

LAMBDA_1=Matrix_N'*UTU*Matrix_N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UTU=d_cheb2(t,M,weight)

    GAMMA=diag(1:M);
    GAMMA=[zeros((M),1) GAMMA];
    
    U(1,1)=0;
        for i=1:M
        U(1,i+1)=i*chebU(t,i-1);%%% size 1\times (M+1)
        end
    
    UTU=U'*U; %%% size (M+1)\times (M+1)
    UTU=weight(t)*UTU;
    
    end

end