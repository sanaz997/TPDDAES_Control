function [flag] = g_fun(phi,phi_en,phi_ex )
if phi>phi_en  && phi<phi_ex
        flag=1;
    else
        flag=0;
end
end

