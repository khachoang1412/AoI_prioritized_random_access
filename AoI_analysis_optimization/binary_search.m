function [xopt] = binary_search(f,x1,x2,TARGET,TOL)
% search the value of x such that f(x) pointwise approaches TARGET from
% below. f(x) is an inceasing function of x.
% xopt = 0;
iter= 30;                       % maximum number of iterations
k1=0;                            % number of iterations
fx=f((x1+x2)/2); %golden_search(f,1e-9,(x1+x2)/2,(x1+x2)/200); % computing values in x points
TARGET = reshape(TARGET,size(fx));
% nDim = length(TARGET);

stop = 0;


while ~stop %(TARGET < fx || fx <TARGET - TOL) && (k1<iter || TARGET < fx)    
    k1=k1+1;
    if k1 > iter
        stop = 1;
    elseif k1 > 1
        fx=f((x1+x2)/2); 
    end
    
    if sum(TARGET > fx) == length(fx)
        x1 = (x1+x2)/2; %set new end of interval        
    else
        x2 = (x1+x2)/2; %replace as new start index
    end
    
    gap = fx - TARGET;
    if sum(gap > 0) == 0 
        if abs(max(gap(fx > 0))) < TOL
            stop = 1;
        end
    end
    
    
end

xopt   = x1;
if f(x1) > 3e4 %sum(TARGET > fx) < length(fx)
    keyboard
end

end
