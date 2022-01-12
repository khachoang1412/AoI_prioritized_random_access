function PLR = PLR_simulation(pClass,x,px,M,G_set,N_rep,iter_max)
% function PLR = computePLR(pClass,x,px,M,G_set,N_rep,iter_max)
% Compute by simulation the PLR of IRSA for a system with K user classes
%
% INPUTS:
%   classDist : fraction of users in each class (set to be 1 if there is one class)
%   x, px     : degree distribution from VN (user) perspective 
%                   (x: degrees, px (K rows): corresponding probability)
%   M         : frame length (#slots)
%   G_set     : set of channel loads for which the PLR is computed
%   N_rep     : number of Monte-Carlo iterations
%   iter_max  : maximum number of decoding iterations
%
% OUTPUTS:
%   PLR       : the PLR for each class and each channel load

x(sum(px) == 0) = [];
px(:,sum(px) == 0) = [];

%% average degree distribution
px_avg = sum(bsxfun(@times,px,pClass(:)),1);

%% Monte-Carlo simulation to compute PLR
PLR = zeros(size(px,1),length(G_set));
for idxG = 1:length(G_set)
    G = G_set(idxG);
    
    transmitted_G = zeros(length(x),1); % count the number of transmitting users per degree
    decoded_G = zeros(length(x),1);     % count the number of decoded users per degree
    
    parfor kkkk = 1:N_rep
        % IRSA encoding and decoding
        [H, L, ptrs] = IRSA_encode(G, M, x, px_avg);
        [~, DECODED] = IRSA_decode(H, L, ptrs, iter_max);
        
        % count the number of transmitting/decoded users per degree
        decoded_tmp = zeros(length(x),1);
        transmitted_tmp = zeros(length(x),1);
        for k = 1:length(x)
            idx = (L == x(k))';
            decoded_tmp(k) = decoded_tmp(k) + sum(DECODED(idx));
            transmitted_tmp(k) =  transmitted_tmp(k) + sum(idx);
        end
        decoded_G =  decoded_G + decoded_tmp;
        transmitted_G =  transmitted_G + transmitted_tmp;
    end
    
    % Compute the PLR
    tmp = 1-decoded_G./transmitted_G;
    tmp(transmitted_G == 0) = 0;
    PLR(:,idxG) = px*tmp;
end