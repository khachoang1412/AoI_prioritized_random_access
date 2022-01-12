function [PLR,G] = PLR_errfloor_multiClass(x, px, M, G, pClass)
% function [PLR,G] = PLR_errfloor_multiClass(x, px, M, G, pClass)
% Compute the approximation of the PLR of multi-class IRSA in the 
% error-floor region according to
%
% [1]  M. Ivanov, F. Brannstrom, A. Graell i Amat, and G. Liva, “Unequal error
% protection in coded slotted ALOHA,” IEEE Wireless Commun. Lett.,
% vol. 5, no. 5, pp. 536–539, Oct. 2016.
%
% INPUTS: 
%   x, px : degree distribution
%   M     : frame length
%   G     : channel load
%   pClass : fraction of users in each class
%
% OUTPUTS: 
%   PLR   : approximate PLR
%   G     : refined channel load

% Number of classes
nClass = length(pClass);    

% Average degree distribution
px_avg = sum(bsxfun(@times,px,pClass(:)),1);

% Strim off the degrees with zero probability
idx = (px_avg==0);
x(idx) = [];
px(:,idx) = [];
px_avg(idx) = [];

%% Approximate the PLR according to [1, Eq.(7)]
PLR = zeros(nClass, length(G));

% PLR per degree for the average degree distribution
[G, PLRtmp, degreeMax] = error_floor_singleClass(x, px_avg, M, G);

% PLR for multiple classes
for idxClass = 1:nClass
    for k = 1:sum(x<=degreeMax)
        if x(k) ==0
            PLR(idxClass,:) = PLR(idxClass,:) + px(idxClass,k);
        else
            PLR(idxClass,:) = PLR(idxClass,:) + PLRtmp(k,:)*px(idxClass,k);
        end
    end
end

end