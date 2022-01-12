function [alfa,beta,z0] = compute_alfa_beta(x,px)
% function [alfa,beta,z0] = compute_alfa_beta(x,px)
% Compute the scaling parameters alpha, beta in finite-length scaling
% analysis according to 
%
% A. Amraoui, A. Montanari, and R. Urbanke, “Analytic determination of
% scaling parameters,” in Proc. IEEE Int. Symp. Inf. Theory (ISIT), Seattle,
% WA, USA, Jul. 2006, pp. 562–566.
%
% INPUTS  : degree distributions x, px
% OUTPUTS : scaling parameters alpha, beta; fixed point z0

rate = 0;

%% Compute distributions
Lambda = [x; px];
Lambda(:,Lambda(2,:) == 0) = [];

Lambda1 = Lambda; 
Lambda1(2,:) = Lambda(2,:).*Lambda(1,:);
Lambda1(1,:) = Lambda(1,:) - 1;
coeff = sum(Lambda1(2,:));
Lambda1(2,:) = Lambda1(2,:)/coeff;

lambda = @(x) degPoly(x,Lambda1);
lambdap = @(x) degPoly_derivative(x,Lambda1);
lambdapp = @(x) degPoly_derivative2(x,Lambda1);

rho = @(x) exp(coeff.*(x-1)./(1-rate));
rhop = @(x) (coeff/(1-rate)).*exp(coeff.*(x-1)./(1-rate));

%% Threshold
thp = thresh_IRSA(x, px);

%% Compute the critical points xfb and yfb 
r1 = @(eps,z) eps.*lambda(z).*(z-1+rho(1-eps.*lambda(z)));
eq = @(z) r1(thp,z);

z0 = linspace(0,1,10000);
eqz0 = eq(z0);
eqoffset = [eqz0(1:end-1) - eqz0(2:end) 0];
eqminus = z0(eqoffset>0);
% eqneg = z0(eqz0<0);
if isempty(eqminus) || sum(eqoffset(z0<eqminus(end)) < 0) == 0  %abs(eq(eqminus(end))) > 1e-4%&& isempty(eqneg)
%    z0 = eps;
   alfa = 0;
   beta = 0;
else
   z0 = eqminus(end);
   
   xfp = thp*lambda(z0);
    yfp = 1- rho(1-xfp);
    xfpb = 1-xfp;
%     yfpb = 1-yfp;

    %% Compute alfa
    Lpoiss = 1/integral(lambda, 0, 1);
    tpa = (rho(xfpb)^2 - rho(xfpb^2) + rhop(xfpb)*(1 - 2*xfp*rho(xfpb)) ...
        - xfpb^2*rhop(xfpb^2))/(Lpoiss*lambda(yfp)^2*rhop(xfpb)^2);
    tpb = thp^2*lambda(yfp)^2 - thp^2*lambda(yfp^2) ...
        - yfp^2*thp^2*lambdap(yfp^2)/(Lpoiss*lambda(yfp)^2);
    alfap = sqrt(tpa + tpb);

    alfa = sqrt(alfap^2 - thp*(1 - thp));

    %% Compute beta
    mmax = 100;
    rhom = @(m) (Lpoiss.^m.*m)./(Lpoiss.*exp(Lpoiss).*factorial(m));
    r = zeros(3,1);
    for i = 2:3
        for m = i:mmax
        for j = i:m
            r(i) = r(i) + (-1)^(i+j)*nchoosek(j-1,i-1)*nchoosek(m-1,j-1)...
                *rhom(m)*(thp*lambda(yfp))^j;
        end
        end
    end
    r2 = r(2); r3 = r(3);

    tbetaa = thp^4*r2^2*(thp*lambdap(yfp)^2*r2 - xfp*(lambdapp(yfp)*r2 + lambdap(yfp)*xfp))^2;
    tbetab = Lpoiss^2*rhop(xfpb)^3*xfp^10*(2*thp*lambdap(yfp)^2*r3 - lambdapp(yfp)*r2*xfp);
    beta = (tbetaa/tbetab)^(1/3);
end
end

%% 
function f = degPoly(y,Lambda)
f = zeros(size(y));
for ii = 1:length(y)
    f(ii) = sum(Lambda(2,:).*y(ii).^(Lambda(1,:)));
end
end

%% 
function f = degPoly_derivative(y,Lambda)
f = zeros(size(y));
for ii = 1:length(y)
    f(ii) = sum(Lambda(2,:).*Lambda(1,:).*y(ii).^(Lambda(1,:)-1).*(Lambda(1,:)>=1));
end
end

%% 
function f = degPoly_derivative2(y,Lambda)
f = zeros(size(y));
for ii = 1:length(y)
    f(ii) = sum(Lambda(2,:).*Lambda(1,:).*(Lambda(1,:)-1).*y(ii).^(Lambda(1,:)-2).*(Lambda(1,:) >= 2));
end
end