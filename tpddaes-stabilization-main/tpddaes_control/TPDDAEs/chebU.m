%%  The file gives in output the value of the Chebyshev polynomial of the second kind with degree n computed in x\in[-1,1]  
function t = chebU( x,n)
    if n==-1
        t = zeros(size(x));
    elseif n == 0                 % use explicit formula when Nch = 0
        t = ones(size(x));
    elseif n == 1             % use explicit formula when Nch = 1
        t = 2*x;
    else                      % use recursive formula when Nch > 1
        ach = 1;
        bch = 2*x;
        for kch = 2 : n
            t = 2 * x .* bch - ach;
            ach = bch;
            bch = t;
        end
    end
end

