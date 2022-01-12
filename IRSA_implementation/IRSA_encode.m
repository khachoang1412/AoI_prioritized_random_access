function [H, L, ptrs] = IRSA_encode(G, M, x, px, pm, type)
% function [H, L, ptrs] = IRSA_encode(G, M, x, px, pm, type)
% Encoding function of IRSA
% 
% INPUT:
%   G     : channel load 
%   M     : frame length (#slots)
%   x, px : degree distribution from VN (user) perspective 
%            (x: degrees, px: corresponding probability)
%   pm    : degree distribution from CN (slot) perspective
%           (if users choose slots uniformly, don't specify pm)
%   type  : distribution of the number of users, either 'fixed' or 'poisson'
%
% OUTPUTS:
%   H     : adjacency matrix, entry (i,j) is 1 if user i transmits in slot j,
%           otherwise 0
%   L     : store the degrees of all users
%   ptrs  : pointers to the slots connected to each user

if ~exist('type', 'var') || isempty(type)
    type = 'poisson';
end

% number of users
if strcmp(type,'fixed')
    U = round(G*M);
elseif strcmp(type,'poisson')
    U = poissrnd(G*M); 
end

% generate the degree of each user
L = randi_distr(x, px, 1, U);
L = sort(L);

% initialize the pointers
maxDegree = max(x);
if ~exist('pn', 'var') || isempty(pm)
    ptrs = zeros(U,maxDegree);
    for k = 1:U
        ptrs(k, :) = randperm(M, maxDegree);
    end   
else
    ptrs = randi_distr(1:M, pm, U, maxDegree);
end
    
% fill in the adjacency matrix and refine the pointers
H = zeros(U, M);
if length(x)>1
    idx0 = 1;
    for k = 1:length(x)-1
        didx = sum(L == x(k)); 
        ptrs(idx0:idx0 - 1 + didx, x(k)+1:end) = NaN;
        idx0 = idx0 + didx;
    end
end
idxs = U*(ptrs-1) + repmat((0:U-1)', 1, max(x)) + 1;
idxs = idxs(isnan(idxs)==0);
H(idxs) = 1;