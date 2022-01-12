function y = randi_distr(x, px, n, m)
% generates an n by m matrix with random entries distributed according to a 
% specified distribution

% px = [0.1 0.1 0.4 0.3 0.1];
% x = [1 2 3 4 5];
tmp = rand(n*m, 1);
cdf = 1./cumsum(px);
idx = (tmp*cdf)<1;
idx(:,2:end) = mod(idx(:,1:end-1) + idx(:,2:end),2);
% use the commented part to generate sorted RVs
% [~, idx] = find(idx);
% y = x(idx);
[idx, ~] = find(idx');
y = x(idx);
if n>1
    y = buffer(y, n);
end



