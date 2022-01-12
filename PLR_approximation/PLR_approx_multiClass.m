function [PLR,PLR_WFEF,PLR_DE,G] = PLR_approx_multiClass(x, px, M, G, U, pClass)
% function [PLR,PLR_WFEF,PLR_DE,G] = PLR_approx_multiClass(x, px, M, G, U, pClass)
% Compute the approximation of the PLR of multi-class IRSA proposed in
%
% [1] K.-H. Ngo, G. Durisi, and A. Graell i Amat, “Age of information in prioritized random access,” in
% 55th Asilomar Conference on Signals, Systems, and Computers, CA, USA, Oct. 2021.
%
% and another one proposed in 
% 
% [2]  A. Munari, “Modern random access: An age of information perspective
% on irregular repetition slotted ALOHA,” IEEE Trans. Commun., vol. 69,
% no. 6, pp. 3572–3585, Jun. 2021.
% 
% and the one obtained by density evolution.
%
% INPUTS: 
%   x, px : degree distribution
%   M     : frame length
%   G     : channel load
%   U     : number of users
%   pClass : fraction of users in each class
%
% OUTPUTS: 
%   PLR   : approximate PLR proposed in [1, Eq.(8)]
%   PLR_WFEF : approximate PLR proposed in [2, Eq.(5)]
%   PLR_DE   : approximate PLR by density evolution
%   G     : refined channel load

nClass = length(pClass);    % number of classes
Gp = G(G>0);

% Average degree distribution
px_avg = sum(bsxfun(@times,px,pClass(:)),1);

% Strim off the degrees with zero probability
x(sum(px) == 0) = [];
px(:,sum(px) == 0) = [];

idx = (px_avg==0);
x(idx) = [];
px(:,idx) = [];
px_avg(idx) = [];

% Decoding threshold for the average degree distribution
if min(x) == 0
    thres = thresh_IRSA(x(x>0),px_avg(x>0)/sum(px_avg(x>0)))/sum(px_avg(x>0));
else
    thres = thresh_IRSA(x,px_avg);
end
    
%% Approximate PLR
PLR = zeros(nClass, length(G));
PLR_WFEF = zeros(nClass, length(G));
PLR_DE = zeros(nClass, length(G));

[~,tmp_DE] = PLR_DE_singleClass(x,px_avg,Gp);
[~, tmp_floor, degreeMax] = PLR_errfloor_singleClass(x, px_avg, M, Gp);
tmp_waterfall = PLR_waterfall(x(x>0),px_avg(x>0)/sum(px_avg(x>0)),M,Gp*sum(px_avg(x>0)),U);
tmp_WFEF = tmp_floor + tmp_waterfall;

for idxClass = 1:nClass
    for k = 1:sum(x<=degreeMax)
        if x(k) == 0 
            PLR(idxClass,G>0) = PLR(idxClass,G>0) + px(idxClass,k);
        elseif x(k) <= 2  %|| ((sum(x == 1) > 0) && sum(px(:,x(x==1)) > 1e-4) > 0)
            PLR(idxClass,G>0) = PLR(idxClass,G>0) + tmp_DE(k,:)*px(idxClass,k);
        else
            PLR(idxClass,G>0) = PLR(idxClass,G>0) + tmp_WFEF(k,:)*px(idxClass,k);
        end
        
        if x(k) == 0 
            PLR_WFEF(idxClass,G>0) = PLR_WFEF(idxClass,G>0) + px(idxClass,k);
            PLR_DE(idxClass,G>0) = PLR_DE(idxClass,G>0) + px(idxClass,k);
        else
            PLR_WFEF(idxClass,G>0) = PLR_WFEF(idxClass,G>0) + tmp_WFEF(k,:)*px(idxClass,k);
            PLR_DE(idxClass,G>0) = PLR_DE(idxClass,G>0) + tmp_DE(k,:)*px(idxClass,k);
        end
    end
end

PLR(:,G>thres+eps) = PLR_DE(:,G>thres+eps); 
end