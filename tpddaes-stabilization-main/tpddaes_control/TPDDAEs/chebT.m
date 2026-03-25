%%  The file gives in output the value of the Chebyshev polynomial of the first kind with degree n computed in x\in[-1,1]  
function t=chebT(x,n)    
    t=cos(n*acos(x));
end