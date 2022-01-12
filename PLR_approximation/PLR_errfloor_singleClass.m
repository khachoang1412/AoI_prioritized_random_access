function [G, PLR, degreeMax] = PLR_errfloor_singleClass(x, px, M, G)
% function [G, PLR, degreeMax] = PLR_errfloor_singleClass(x, px, M, G)
% Compute the approximation of the PLR of single-class IRSA in the 
% error-floor region according to
%
% [1] M. Ivanov, F. Brannstrom, A. Graell i Amat, and P. Popovski, “Broadcast
% coded slotted ALOHA: A finite frame length analysis,” IEEE Trans.
% Commun., vol. 65, no. 2, pp. 651–662, Feb. 2017.
%
% INPUTS: 
%   x, px : degree distribution
%   M     : frame length
%   G     : channel load
%
% OUTPUTS: 
%   G     : refined channel load
%   PLR   : approximate PLR, store the per-degree PLR concatenated by the
%               overall PLR
%   degreeMax : maximum degree

%% If there are only degrees 0 and 1 (slotted ALOHA), PLR is known in closed form
if max(x) == 1
    degreeMax = 1;
    if sum(x==0) == 0
        plr_deg1 = 1 - (1-1/M).^(M.*G-1);
        PLR = [plr_deg1; plr_deg1];
    else
        plr_deg1 = 1 - (1-1/M).^(M.*G*px(x==1)-1);
        PLR = [ones(size(plr_deg1)); plr_deg1; px(x==0) + px(x==1)*plr_deg1];
    end
else

%% Sort the degrees in ascending order
[x, idx] = sort(x);
px = px(idx);

%% Strim off the degrees with zero probability
idx = (px==0);
x(idx) = [];
px(idx) = [];

%% Load the minimal stopping sets and initialize some parameters
ss_list31       
% ss_list142 % bigger set

degreeMax = size(v,2);
px_full = zeros(1, degreeMax);
for k = 1:degreeMax
    idx = find(x == k);
    if isempty(idx) ==0
        px_full(k) = px(idx);
    end
end
U = round(M*G); % number of users
U = unique(U);  
G = U/M;

%% Approximation according to [1, Eq.(30)]
ss_num = length(c);   % number of sets of stopping sets  

% Compute nu for each set of stopping sets
nu = sum(v,2)';

% Compute alpha and beta 
alpha = ones(size(c));
beta = zeros(size(c));

for k = 1:ss_num
    alpha(k) = factorial(nu(k));    
    beta(k)= c(k)*nchoosek(M, mus(k));
    for kk = 1:degreeMax
        alpha(k) = alpha(k)*px_full(kk)^v(k, kk)/factorial(v(k, kk));        
        beta(k) = beta(k)/(nchoosek(M, kk)^v(k, kk));
    end
end

% Compute gamma
gamma = zeros(ss_num, length(U));
for k = 1:ss_num
    for kk = 1:length(U)
%         if U(kk) < nu(k)
%             keyboard
%         end
        gamma(k, kk) = alpha(k)*beta(k)*nchoosek(U(kk), nu(k));
    end
end

% Put them together to compute the PLR
PLR = zeros(sum(x<=degreeMax)+1, length(U));

for k = 1:sum(x<=degreeMax)
    if x(k) ==0
        PLR(k,:) = 1;
        PLR(end,:) = PLR(end,:) + px(k);
    else
        for kk = 1:ss_num
            PLR(k, :) = PLR(k, :) + v(kk, x(k))*gamma(kk, :);
        end
        PLR(k, :) = PLR(k,:)./U/px(k);
        PLR(end,:) = PLR(end,:) + PLR(k,:)*px(k);
    end
end

end