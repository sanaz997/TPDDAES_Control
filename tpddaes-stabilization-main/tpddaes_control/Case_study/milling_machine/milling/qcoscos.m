function [ qval ] = qcoscos(t,omega,N,pp,q,phi_en,phi_ex)
qval=0;
if phi_angle(1,t,omega,N)==0    
else
    qval=qval+pp*g_fun(phi_angle( 1,t,omega,N ),phi_en,phi_ex )*cos(phi_angle(1,t,omega,N))^2*sin(phi_angle(1,t,omega,N))^(q-1);
end

if phi_angle(2,t,omega,N)==0    
else
    qval=qval+pp*g_fun(phi_angle(2,t,omega,N ),phi_en,phi_ex )*cos(phi_angle(2,t,omega,N))^2*sin(phi_angle(2,t,omega,N))^(q-1);
end


end

