% Evaluate the packet loss rate (PLR) and age-violation probability (AVP) 
% for a system with 2 classes of users operating according to IRSA
%
% [1] K.-H. Ngo, G. Durisi, and A. Graell i Amat, “Age of information in prioritized random access,” in
% 55th Asilomar Conference on Signals, Systems, and Computers, CA, USA, Oct. 2021.

clear 
close all

addpath(genpath('..'));

%% PARAMETERS:
U = 4000;                   % number of users
nClass = 2;                 % number of classes
pClass = [0.2 0.8];         % fraction of users in each class
peakAge = [7.5e4 4.5e4];    % AoI thresholds

% Degree distribution
x = [1 2 3];
px = [0 .5 .5;
      0 .5 .5];
x(sum(px) == 0) = [];
px(:,sum(px) == 0) = [];
px = bsxfun(@rdivide,px,sum(px,2));

% Average degree distribution and the decoding threshold
px_avg = pClass*px;
thres = thresh_IRSA(x, px_avg);
P0 = x*px.';

% Maximum number of decoding iterations in IRSA
iter_max = 40;

% Sets of considered frame lengths and channel loads
M_set = 100; %:200:2000;    % frame length 
G_set = [0.05 [.1:.1:3.5]/((pClass(:).'*P0(:)))];    % channel load    
G_set = G_set(G_set >= 0.05);

%% IRSA, based on [1]
pktReadyProb_set = zeros(length(G_set),length(M_set));
ProbActive_set = zeros(length(G_set),length(M_set));
Power = zeros(length(G_set),length(M_set)); % denoted by \Phi in [1]

PLR_approx = zeros(nClass,length(G_set),length(M_set));
PLR_WFEF = zeros(nClass,length(G_set),length(M_set));
PLR_DE = zeros(nClass,length(G_set),length(M_set));
PLR = zeros(nClass,length(G_set),length(M_set));

AVP = zeros(nClass,length(G_set),length(M_set));
AVP_approx = zeros(nClass,length(G_set),length(M_set));
AVP_WFEF = zeros(nClass,length(G_set),length(M_set));
AVP_DE = zeros(nClass,length(G_set),length(M_set));
    
for idxG = 1:length(G_set)    
G = G_set(idxG)

% set the number of Monte-Carlo iterations
if G < 0.5         % fewer iterations if the channel load is low
    N_rep = 1e3;
else
    N_rep = 1e4;
end

for idxM = 1:length(M_set)  
    M = M_set(idxM);
    
    % Prob that a user has an update in a frame, denoted by \sigma in [1]
    pktReadyProb_set(idxG,idxM) = M*G/U;
    pktReadyProb = pktReadyProb_set(idxG,idxM); 
    
    % Prob that a user has an update in a slot, denoted by \mu in [1]
    ProbActive_set(idxG,idxM) = 1-(1-pktReadyProb)^(1/M);
    ProbActive = ProbActive_set(idxG,idxM); 
            
    % Avg number of transmitted packets, denoted by \Phi in [1]
    Power(idxG,idxM) = pktReadyProb_set(idxG,idxM)*U/M*(pClass(:).'*P0(:));

    % Exact PLR
    PLR(:,idxG,idxM) = PLR_simulation_2(pClass,x,px,M,U,pktReadyProb,N_rep,iter_max);
    
    % Approximate PLR
    [plrApprox,plrWFEF,plrDE] = PLR_approx_multiClass(x,px,M,G,U,pClass);
    PLR_approx(:,idxG,idxM) = min(plrApprox,1);     % [1, Eq.(8)]
    PLR_WFEF(:,idxG,idxM) = min(plrWFEF,1);         % [1, Eq.(7)], according to Munari
    PLR_DE(:,idxG,idxM) = min(plrDE,1);             % by density evolution

    % Age violation probability
    for idxClass = 1:nClass
        AVP(idxClass,idxG,idxM) = computeAVP(M,...
            ProbActive,peakAge(idxClass),PLR(idxClass,idxG,idxM));
        AVP_approx(idxClass,idxG,idxM) = computeAVP(M,...
            ProbActive,peakAge(idxClass),PLR_approx(idxClass,idxG,idxM));
        AVP_WFEF(idxClass,idxG,idxM) = computeAVP(M,...
            ProbActive,peakAge(idxClass),PLR_WFEF(idxClass,idxG,idxM));
        AVP_DE(idxClass,idxG,idxM) = computeAVP(M,...
            ProbActive,peakAge(idxClass),PLR_DE(idxClass,idxG,idxM));
    end
end 
end

%% Plot figure

% PLR vs Power (Fig 2(a), 2(b) in [1])
figure;
semilogy(Power,PLR,'r-')
hold on
semilogy(Power,PLR_WFEF,'b--')
semilogy(Power,PLR_DE,'k:')
semilogy(Power,PLR_approx,'m-.')
ylim([1e-5 1])
grid on
xlabel('$\Phi$ (packets/slot)','interpreter','latex')
ylabel('PLR','interpreter','latex')
legend('Simulation','$P_\ell \approx P_{\rm WF} + P_{\ell, {\rm EF}}$',...
    'P_\ell \approx P_{\ell,{\rm DE}}','Proposed approx.')

% AVP vs Power (Fig 2(c), 2(d) in [1])
figure
semilogy(Power,AVP,'r-')
hold on
semilogy(Power,AVP_WFEF,'b--')
semilogy(Power,AVP_DE,'k:')
semilogy(Power,AVP_approx,'m-.')
grid on
ylim([1e-6 1])
xlabel('$\Phi$ (packets/slot)','interpreter','latex')
ylabel('AVP','interpreter','latex')

%% Slotted ALOHA (see Proposition 4 in Munari (2021)
SA_AVP = zeros(nClass,length(G_set),length(M_set));
SA_peak_tolerable_AoI = zeros(nClass,length(G_set),length(M_set));
SA_Power = zeros(length(G_set),length(M_set));

for idxG = 1:length(G_set)    
G = G_set(idxG);
for idxM = 1:length(M_set) 
    M = M_set(idxM);
    pktReadyProb = M*G/U; 
    ProbActive = 1-(1-pktReadyProb)^(1/M); 

    SA_Power(idxG,idxM) = U*ProbActive;
    SA_PLR = ProbActive*(1-ProbActive)^(U-1);
    SA_AVP(:,idxG,idxM) = (1-SA_PLR).^(peakAge-1);
end 
end