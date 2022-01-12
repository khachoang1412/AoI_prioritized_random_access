function [PLR,FER] = PLR_waterfall(x,px,M,G,U)
% function [PLR,FER] = PLR_waterfall(x,px,M,G,U)
% Compute the approximation of the PLR and FER (frame error rate) of IRSA
% in the waterfall region according to
%
% [1] A. Graell i Amat and G. Liva, “Finite-length analysis of irregular
% repetition slotted ALOHA in the waterfall region,” IEEE Commun. Lett.,
% vol. 22, no. 5, pp. 886–889, May 2018.
%
% INPUTS: 
%   x, px : degree distribution
%   M     : frame length
%   G     : channel load
%   U     : number of users
%
% OUTPUTS: packet loss rate PLR, frame error rate FER

%% Compute the scaling parameters alfa and beta
[alfa,beta] = compute_alfa_beta(x,px);

%% Compute the decoding threshold
thp = thresh_IRSA(x, px);

%% Compute the asymptotic PLR when G -> 1
iter = 1e3;
coeff = sum(px.*x);
lambda_l = px.*x/coeff;
pk = 1;
for kk = 1:iter
   qk = sum(lambda_l(1:end).*pk.^(x(1:end)-1));
   pk = 1 - exp(-1*coeff*qk); % here G = 1
end
PLR_G1 = sum(px.*pk.^x);

%% Compute approx. FER and PLR 
rate = M.*G/U;

% [1, eq. (10)]
FER = qfunc(sqrt(M).*(thp - beta*M^(-2/3)-G)./sqrt(alfa^2+G.*(1-rate))); 

% [1, eq. (12)]
PLR = PLR_G1*FER;
end