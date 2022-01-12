function [PLR,PLR_per_degree] = PLR_DE_singleClass(x,px,G)
% function [PLR,PLR_per_degree] = PLR_DE(x,px,G)
% Compute the asymptotic PLR when the framelength goes to infinity based on
% density evolution (DE), according to Proposition 1 in 
% 
% M. Ivanov, F. Brannstrom, A. Graell i Amat, and G. Liva, “Unequal error
% protection in coded slotted ALOHA,” IEEE Wireless Commun. Lett.,
% vol. 5, no. 5, pp. 536–539, Oct. 2016.
%
% INPUTS: 
%   x, px : degree distribution (there can be degree 0)
%   G     : channel load
%
% OUTPUTS:
%   PLR   : overall PLR
%   PLR_per_degree : PLR for each degree

%% Regularize the degree distribution in case there is degree 0
x0 = x;
px0 = px;

if ~isempty(px0(x0==0))
    offset = px0(x0==0);
else 
    offset = 0;
end

probActive = sum(px(x0>0));  % probability for a user to be active
px(x0==0) = [];
px = px/probActive;

x(x==0) = [];

G = G*probActive;

%% Density evolution
iter = 1e4;
PLR = zeros(size(G));
eta = zeros(size(G));   % probability that an edge remains unknown  

flag = 1;
k = 1;

coeff = sum(px.*x);     % Lambda(1)
lambda = px.*x/coeff;

% Overall PLR
while (flag)&&(k<=length(G))
    pk = 1;
    for kk = 1:iter
       qk = sum(lambda(1:end).*pk.^(x(1:end)-1));
       pk = 1 - exp(-G(k)*coeff*qk);
       zk = sum(px.*pk.^x);
    end
    PLR(k) = zk;
    eta(k) = pk;
    k = k + 1;
end
PLR = offset + PLR*probActive;

% Per-degree PLR
PLR_per_degree = zeros(length(x0), length(G));
for kk = 1:length(x0)
    if x0(kk) == 0
        PLR_per_degree(kk, :) = 1;
    else
        PLR_per_degree(kk, :) = (eta).^x0(kk);
    end
end
end