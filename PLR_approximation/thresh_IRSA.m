function Gmax = thresh_IRSA(x, px)
% function Gmax = thresh_IRSA(x, px)
% Compute the decoding threshold for IRSA with degree distribution x, px

if px(x == 0) > 0 
    Gmax = 0;
else
    coeff = sum(px.*x);         % derivative of the degree distribution evaluated at 1
    lambda_l = px.*x/coeff;     % probability that an edge is incident to a degree-l VN

    f1 = @(y) sum(lambda_l.*y.^(x-1)); % edge perspective VN degree polynomial
    f2 = @(y, G) exp( -G*coeff*(1-y)); % exp(G*(derivative at y))

    dp = 1e-3;
    p = 0:dp:1;             % search for p in (0,1]
    qk = zeros(size(p));    % edge perspective degree polynomial evaluated at p
    for k = 1:length(p)
       qk(k) =  f1(p(k));
    end

    f3 = @(G) sum(p - 1 + f2((1 - qk),G) < 0)-1;
    Gmax = fzero(f3, 1);    % find the largest G such that p < 1 - f2((1 - qk),G)
end
