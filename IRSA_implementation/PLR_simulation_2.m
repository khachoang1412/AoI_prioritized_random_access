function PLR = PLR_simulation_2(pClass,x,px,M,U,pktReadyProb,N_rep,iter_max)
% function PLR = PLR_simulation_2(pClass,x,px,M,U,pktReadyProb,N_rep,iter_max)
% Compute by simulation the PLR of IRSA for a system with K user classes
%
% Different from PLR_simulation(), here we allow for degree 0.
%
% INPUTS:
%   classDist : fraction of users in each class (set to be 1 if there is one class)
%   x, px     : degree distribution from VN (user) perspective 
%                   (x: degrees, px (K rows): corresponding probability)
%   M         : frame length (#slots)
%   U         : number of users 
%   pktReadyProb : probability that a user has a new packet in a frame
%   N_rep     : number of Monte-Carlo iterations
%   iter_max  : maximum number of decoding iterations
%
% OUTPUTS:
%   PLR       : the PLR for each class

if numel(pktReadyProb) == 1
    pktReadyProb = pktReadyProb*ones(length(pClass),1);
end

x(sum(px) == 0) = [];
px(:,sum(px) == 0) = [];

% check if there is degree 0
if min(x) == 0
    offset = px(:,x == 0);
    px(:,x == 0) = [];
    x(x == 0) = [];
    px = bsxfun(@rdivide,px,sum(px,2));
    G_set = (U.*pClass(:)'*bsxfun(@times,pktReadyProb,1-offset)/M);
else
    offset = 0;
    G_set = (U.*pClass(:)'*pktReadyProb/M);
end

% average degree distribution
px_avg = sum(bsxfun(@times,px,pClass(:)),1);

% Monte-Carlo simulation to compute PLR
PLR = zeros(size(px,1),size(pktReadyProb,2));

for idxG = 1:length(G_set)
    G = G_set(idxG);
    
    decoded_G = zeros(length(x),1);
    transmitted_G = zeros(length(x),1);
    
    parfor kkkk = 1:N_rep
        % IRSA encoding and decoding
        [H, L, ptrs] = IRSA_encode(G, M, x, px_avg,[],'poisson');
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
    if sum(isnan(PLR(:))) > 0
        keyboard
    end
end

PLR = bsxfun(@plus,bsxfun(@times,PLR,1-offset),offset);