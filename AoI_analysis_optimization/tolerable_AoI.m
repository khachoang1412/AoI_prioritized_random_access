function tolerable_AoI = tolerable_AoI(M,pktReadyProb,AVP,PLR)
% function tolerable_AoI = tolerable_AoI(M,pktReadyProb,AVP,PLR)
% Compute the AoI threshold which the AoI value does not exceed with 
% certain probability 
%
% INPUTS: 
%   M 	: frame length
%   pktReadyProb : prob that a user has a new update in a frame
%   AVP : target age violation probability
%   PLR : packet loss rate
%
% OUTPUTS: the AoI threshold for each class

nClass = length(AVP);        % number of classes

tolerable_AoI = zeros(nClass,1);    

ProbActive = (1-(1-pktReadyProb).^(1/M));   % prob. that a user is active in a slot (\mu)

%% Compute tolerable AoI
% function to compute the AVP given age value (y) and packet loss rate (plr)
beta = @(y) round((y-2*M)./M);
alpha = @(y) mod((y-2*M),M);

% compute the age value such that the AVP does not exceed the target for
% each class 
for idxClass = 1:nClass
    if ProbActive(idxClass) > 0
        if AVP(idxClass) == 1
            tolerable_AoI(idxClass) = 0;
        elseif AVP(idxClass) == 0
            tolerable_AoI(idxClass) = inf;
        else 
            a = @(y) 1-((1-ProbActive(idxClass)).^(1+alpha(y)) - ...
                (1-ProbActive(idxClass))^M)./(1 - (1-ProbActive(idxClass))^M);
            f1 = @(y,plr) (1-(1-plr).*pktReadyProb(idxClass)).^beta(y)...
                .*(1-a(y).*(1-plr).*pktReadyProb(idxClass));

            f2 = @(plr) fzero(@(y) f1(y,plr)-AVP(idxClass),...
                2*M-1+log(AVP(idxClass))./log(1-ProbActive(idxClass)));
            tolerable_AoI(idxClass) = f2(PLR(idxClass));
        end
    else
        tolerable_AoI(idxClass) = inf;
    end
end