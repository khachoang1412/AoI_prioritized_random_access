function PLR = PLR_DE_multiClass(x, px, G, pClass)
% function PLR = PLR_DE_multiClass(x, px, G, pClass)
% Compute the asymptotic PLR when the framelength goes to infinity for 
% multi-class IRSA based on density evolution (DE)
%
% INPUTS: 
%   x, px : degree distribution (there can be degree 0)
%   G     : channel load
%   pClass : fraction of users in each class
%
% OUTPUTS:
%   PLR   : approximate PLR for each class

% number of classes
nClass = length(pClass);    

% average degree distribution
px_avg = sum(bsxfun(@times,px,pClass(:)),1);

% Strim off the degrees with zero probability
idx = (px_avg==0);
x(idx) = [];
px(:,idx) = [];
px_avg(idx) = [];

% Approximate PLR
PLR = zeros(nClass, length(G));
[~,tmp_DE] = PLR_DE_singleClass(x,px_avg,G(G>0));
    
for idxClass = 1:nClass
for k = 1:sum(x<=degreeMax)
    if x(k) == 0
        PLR(idxClass,G>0) = PLR(idxClass,G>0) + px(idxClass,k);
    else
        PLR(idxClass,G>0) = PLR(idxClass,G>0) + tmp_DE(k,:)*px(idxClass,k);
    end
end
end

end