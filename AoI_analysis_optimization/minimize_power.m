function data = minimize_power(M,U,x,pClass,peakAge,AVP,PLR_mode)
% function data = minimize_power(M,U,x,pClass,peakAge,AVP,PLR_mode)
% Optimize the degree distribution and the activation probability for two 
% classes of users following IRSA in order to minimize the average number 
% of transmitted packets per slot while guaranteeing that the AoI does not 
% exceed a given threshold with a certain probability.
% 
% [1] K.-H. Ngo, G. Durisi, and A. Graell i Amat, “Age of information in prioritized random access,” in
% 55th Asilomar Conference on Signals, Systems, and Computers, CA, USA, Oct. 2021.
%
% INPUTS: 
%   M       : frame length
%   U       : number of users
%   x       : degrees
%   pClass  : fraction of users in each class
%   peackAge : age threshold
%   AVP     : target age violation probability
%   PLRmode : how to compute the PLR? 'simulation' or 'approximation'?
%
% OUTPUTS:
%   data    : store the inputs and the optimization results

addpath(genpath('..'));

tic

% DEBUG mode?
DEBUG = 1;

if DEBUG
    M = 100;    
    U = 4000;  
    pClass = [0.2 0.8];
    x = [0 2 3];      
    AVP = [1e-3 1e-3];
    peakAge = [7.5e4 4.5e4];
    PLR_mode = 'approximation';
    ini_mode = 'random';
    timeout = 2*60;
end

warning off

%% PARAMETERS
% number of classes (currently we allow for 2 classes only)
nClass = length(pClass);

% precompute a matrix that will be repeatedly used
A = kron(eye(nClass),ones(1,length(x)-1));
A = [zeros(size(A,1),1) A];

% number of Monte-Carlo iterations for IRSA
N_rep = 1e4;

% number of IRSA decoding iterations
iter_max = 40;

%% INITIALIZATION
% multiple initializations for the degree distributions
% two options: 1) sampling the search space with a step
%              2) sampling the search space randomly

ini_mode = 'random'; % 'step' or 'random'

% Note: if you choose 'step', modify the set of initialization below 
%   according to the number of degrees
step = 0.1;
iini = 1;
for lambda11 = 0:step:1
for lambda12 = 0:step:min(1-lambda11,1)
% for lambda13 = 0:step:1-lambda11-lambda12
for lambda21 = 0:step:1
for lambda22 = 0:step:min(1-lambda21,1)
% for lambda23 = 0:step:1-lambda21-lambda22
    px_ini(:,:,iini) = [lambda11    lambda12; 
           lambda21    lambda22];
    iini = iini + 1;
end
end
end
end

% If you choose 'random', set a timeout 
timeout = 3600; 

%% OPTIMIZATION
G_opt = 0;
pktReadyProb_opt = 0;
mu_opt = 0;
px_opt = zeros(nClass,length(x),1);
power_opt = 0;
PLR_opt = zeros(nClass,1);
peakAge_opt = zeros(nClass,1);
AVP_opt = zeros(nClass,1);

ii = 1;
idxIni = 1;

flag = 1;
while flag
    % initialize the degree distribution
    if strcmp(ini_mode,'step')
        px0 = px_ini(:,:,idxIni)';
    elseif strcmp(ini_mode,'random')
        px0 = randfixedsum(length(x),nClass,1,0,1).';
        
        % or you can play with the below to refine the range of sampling
%         px0 = zeros(nClass,length(x));
%         for idxX = 1:length(x)
%             if x(idxX) >= 3
%                 px0(:,idxX) = rand(nClass,1);
%             else
%                 px0(:,idxX) = rand(nClass,1)*0.5;
%             end
%         end
%         px0 = bsxfun(@rdivide,px0,sum(px0,2));

    end
    px0 = px0(:,1:end-1)';
    
    % how to compute the PLR?
    if strcmp(PLR_mode,'approximation')
        PLRcomputation = @(pktReadyProb,px) PLR_approx_multiClass(x,[reshape(px,length(x)-1,[]).' ...
                    1-sum(reshape(px,length(x)-1,[]).',2)],M,pktReadyProb*U/M,U,pClass);
    elseif strcmp(PLR_mode,'simulation')
        PLRcomputation = @(pktReadyProb,px) PLR_simulation(pClass,x,[reshape(px,length(x)-1,[]).' ...
                    1-sum(reshape(px,length(x)-1,[]).',2)],M,...
                    pktReadyProb*U/M,N_rep,iter_max);
    end
    
    % function to compute the tolerable AoI given target AVP
    f = @(pktReadyProb,px) tolerable_AoI(M,pktReadyProb*ones(1,nClass),...
            AVP,PLRcomputation(pktReadyProb,px));
    
    % set up the AoI constraint: the tolerable AoI must be lower than
    % the target threshold
    nonlcon = @(optVar) f(optVar(1),optVar(2:end)) - peakAge(:);
    
    % innitialize the activation probability \mu such that the AoI
    % constraint is feasible
    px0 = [0.0164; px0(:)];
    check = 1;
    if sum(nonlcon(px0) > 0) > 0
        check = 0;
        for gtmp = .3:0.02:.9
            px0(1) = gtmp*M/U;
            if sum(nonlcon(px0) > 0) == 0
                check = 1;
                break;
            end
        end
    end
    
    if check    % feasible :), go ahead!
        % function compute the power (i.e., avg number of transmitted
        % packets per slot
        power = @(optVar) (U/M)*optVar(1)*((x*([reshape(optVar(2:end),length(x)-1,[]).' ...
                   1-sum(reshape(optVar(2:end),length(x)-1,[]).',2)]).'))*pClass(:);
        
        % upper limit for the optimization variables
        upperLimit = ones(size(px0));
        upperLimit(1) = M/U;
        
        % you can play around with the limits for degree distribution
%         upperLimit = ones(nClass,length(x)-1);
%         upperLimit(:,x(1:end-1) < 3) = 0.4 + rand*0.1;
%         upperLimit = upperLimit.';
%         upperLimit = [M/U; upperLimit(:)];
        
        % optimization
        [opt_result, obj_opt,~,output] = fminsearchcon(power,px0,zeros(size(px0)),upperLimit,...
                A,ones(nClass,1),nonlcon);
            
        % collect the results
        if ~output.infeasible
            pktReadyProb_opt(ii) = opt_result(1);
            G_opt(ii) = pktReadyProb_opt(ii)*U/M;
            mu_opt(ii) = 1-(1-pktReadyProb_opt(ii)).^(1/M);
            
            opt_result = reshape(opt_result(2:end),length(x)-1,[]).';
            px_opt(:,:,ii) = [opt_result 1-sum(opt_result,2)];
            
            power_opt(ii) = obj_opt; 
            plr = PLR_simulation(pClass,x,px_opt(:,:,ii),M,pktReadyProb_opt(ii)*U/M,1e4,iter_max);
            PLR_opt(:,ii) = plr(:);
            peakAge_opt(:,ii) = tolerable_AoI(M,pktReadyProb_opt(ii)*ones(1,nClass),AVP,plr);
            for idxClass = 1:nClass
                AVP_opt(idxClass,ii) = computeAVP(M,mu_opt(ii),...
                    peakAge(idxClass),plr(idxClass));
            end
            ii = ii + 1;
        end
    end
    
    idxIni = idxIni + 1;
    
    % check termination condition
    if (strcmp(ini_mode,'step') && idxIni > size(px_ini,3)) ...
            || (strcmp(ini_mode,'random') && toc > timeout)
        flag = 0;
    end
end

%% Save the results
data.M = M;
data.U = U;
data.x = x;
data.nClass = nClass;
data.pClass = pClass;
data.AVP = AVP;
data.peakAge = peakAge;


[power_opt, ia] = min(power_opt);

data.power_opt = power_opt;
data.pktReadyProb_opt = pktReadyProb_opt(ia);
data.G_opt = G_opt(ia);
data.mu_opt = mu_opt(ia);
data.px_opt = px_opt(:,:,ia);
data.PLR_opt = PLR_opt(:,ia);
data.peakAge_opt = peakAge_opt(:,ia);
data.AVP_opt = AVP_opt(:,ia);

data.simtime = toc;

if DEBUG
    keyboard
else
    if ii > 1
        filename = ['Opt_m_' num2str(M) '_n_' num2str(U) '_pClass1_' num2str(pClass(1))...
            '_AVP_' num2str(AVP(1)) '_' num2str(AVP(end)) ...
            '_xMax_' num2str(max(x)) '_power_' num2str(power_opt) '.mat'];

        save(filename, 'data', '-v7.3');
    else
        disp('infeasible');
    end 
end
end