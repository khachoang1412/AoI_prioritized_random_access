function AVP = computeAVP(M,probActive,peakAge,PLR)
% function AVP = computeAVP(M,probActive,peakAge,PLR)
% Compute the age violation probability (AVP) according to Proposition 2 in
%
% [1] K.-H. Ngo, G. Durisi, and A. Graell i Amat, “Age of information in prioritized random access,” in
% 55th Asilomar Conference on Signals, Systems, and Computers, CA, USA, Oct. 2021.
%
% INPUTS: 
%   M           : frame length
%   probActive  : prob that a user has a new update in a slot, denoted by
%                   \mu in [1]
%   peakAge     : AoI threshold
%
% OUTPUT: 
%   AVP         : age violation probability

AVP = zeros(length(peakAge),length(PLR));

if peakAge <= 2*M
    warning('Impossible to achieve this peak age!')
else 
    warning off;
     
    % compute the AVP according to [1, Prop. 2]
    for iii = 1:length(peakAge)
        beta = round((peakAge(iii)-2*M)./M);
        alpha = mod((peakAge(iii)-2*M),M);
        a = 1-((1-probActive).^(1+alpha) - (1-probActive)^M)./(1 - (1-probActive)^M);
        x = (1 - (1-probActive)^M).*(1-PLR);

        AVP(iii,:) = (1-x).^beta.*(1-a*x);
    end
    AVP(AVP<0) = 0;
    AVP(AVP>1) = 1;
end