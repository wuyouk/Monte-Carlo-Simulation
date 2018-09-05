% ============================================================
%    Calculating Greeks via Simulation (Black-Merton-Scholes)
%
%           (1) w/o Common Random Numbers
%           (2) w/ Common Random Numbers
%           (3) via Pathwise Estimator
%           (4) via Likelihood Ratio
%          
% -------------------------------------------------------------

%% Set up parameters
s0 = 100;
K = 115;
sig = 0.28;
T = 1;
r = 0.04;
q = 0.015;
H = 130;
delS = 0.1;
delSig = 0.001;
nN = 500000;
m = 12;
dt = T / m;

%% Simulation
z1 = randn(nN,1);
z2 = randn(nN,1);
z3 = randn(nN,1);

% ===================================
% ===================================
%
%         Delta & Gamma
%
% -----------------------------------
tmp1U = (s0+delS)*exp((r-q-sig^2/2)*T);
tmp1  =  s0      *exp((r-q-sig^2/2)*T);
tmp1D = (s0-delS)*exp((r-q-sig^2/2)*T);
tmp2 = sig*sqrt(T);

% ---------------------------------
% (1) W/ Common Random Numbers
s = tmp1*exp(tmp2*z1);
su = tmp1U*exp(tmp2*z1);
sd = tmp1D*exp(tmp2*z1);

% (a) Digital Call
tmp = exp(-r*T) * (s > K);
tmpU = exp(-r*T) * (su > K);
tmpD = exp(-r*T) * (sd > K);

%delDC_CRN1 = mean((tmpU-tmp))/delS;
delDC_CRN = mean((tmpU-tmpD))/(2*delS);
gamDC_CRN  = mean((tmpU-2*tmp+tmpD))/(delS^2);

% (b) Up and Out Call
tmp = exp(-r*T) * (s < H) .* max(s - K, 0);
tmpU = exp(-r*T) * (su < H) .* max(su - K, 0);
tmpD = exp(-r*T) * (sd < H) .* max(sd - K, 0);

%delUOC_CRN1 = mean((tmpU-tmp))/delS;
delUOC_CRN = mean((tmpU-tmpD))/(2*delS);
gamUOC_CRN  = mean((tmpU-2*tmp+tmpD))/(delS^2);

% (c) Average Option Call
S_aoc = zeros(nN,m);
Su_aoc = zeros(nN,m);
Sd_aoc = zeros(nN,m);
z_aoc = randn(nN, m);
zu_aoc = randn(nN, m);
zd_aoc = randn(nN, m);
for i = 1:nN
    S_previous = s0;
    Su_previous = s0+delS;
    Sd_previous = s0-delS;
    for ite = 1:m
        if ite > 1
            S_previous = S_aoc(i, ite - 1);
            Su_previous = Su_aoc(i, ite - 1);
            Sd_previous = Sd_aoc(i, ite - 1);
        end
        S_aoc(i, ite) = S_previous * exp((r - q - sig * sig / 2) * dt + sig * sqrt(dt) * z_aoc(i, ite));
        Su_aoc(i, ite) = Su_previous * exp((r - q - sig * sig / 2) * dt + sig * sqrt(dt) * z_aoc(i, ite));
        Sd_aoc(i, ite) = Sd_previous * exp((r - q - sig * sig / 2) * dt + sig * sqrt(dt) * z_aoc(i, ite));
    end
end


tmp = exp(-r*T) * max(mean(S_aoc, 2) - K, 0);
tmpU = exp(-r*T) * max(mean(Su_aoc, 2) - K, 0);
tmpD = exp(-r*T) * max(mean(Sd_aoc, 2) - K, 0);
delAOC_CRN = mean((tmpU-tmpD))/(2*delS);
gamAOC_CRN  = mean((tmpU-2*tmp+tmpD))/(delS^2);

% ----------------------------------
% (2) W/O Common Random Numbers
s =  tmp1*exp(tmp2*z1);
su = tmp1U*exp(tmp2*z2);
sd = tmp1D*exp(tmp2*z3);

% (a) Digital Call
tmp = exp(-r*T) * (s > K);
tmpU = exp(-r*T) * (su > K);
tmpD = exp(-r*T) * (sd > K);

%delDC_tilde1 = mean((tmpU-tmp ))/delS;
delDC_tilde = mean((tmpU-tmpD))/(2*delS);
gamDC_tilde  = mean((tmpU-2*tmp+tmpD))/(delS^2);

% (b) Up and Out Call
tmp = exp(-r*T) * (s < H) .* max(s - K, 0);
tmpU = exp(-r*T) * (su < H) .* max(su - K, 0);
tmpD = exp(-r*T) * (sd < H) .* max(sd - K, 0);

%delUOC_tilde1 = mean((tmpU-tmp ))/delS;
delUOC_tilde = mean((tmpU-tmpD))/(2*delS);
gamUOC_tilde  = mean((tmpU-2*tmp+tmpD))/(delS^2);

% (c) Average Option Call
Su_aoc = zeros(nN,m);
Sd_aoc = zeros(nN,m);
for i = 1:nN
    Su_previous = s0+delS;
    Sd_previous = s0-delS;
    for ite = 1:m
        if ite > 1
            Su_previous = Su_aoc(i, ite - 1);
            Sd_previous = Sd_aoc(i, ite - 1);
        end
        Su_aoc(i, ite) = Su_previous * exp((r - q - sig * sig / 2) * dt + sig * sqrt(dt) * zu_aoc(i, ite));
        Sd_aoc(i, ite) = Sd_previous * exp((r - q - sig * sig / 2) * dt + sig * sqrt(dt) * zd_aoc(i, ite));
    end
end

tmp = exp(-r*T) * max(mean(S_aoc, 2) - K, 0);
tmpU = exp(-r*T) * max(mean(Su_aoc, 2) - K, 0);
tmpD = exp(-r*T) * max(mean(Sd_aoc, 2) - K, 0);
delAOC_tilde = mean((tmpU-tmpD))/(2*delS);
gamAOC_tilde  = mean((tmpU-2*tmp+tmpD))/(delS^2);


% ---------------------------------
% (3) Pathwise Estimator

% (a) Digital Call
delDC_pathwise = NaN;
gamDC_pathwise = NaN;

% (b) Up and Out Call
delUOC_pathwise = NaN;
gamUOC_pathwise = NaN;

% (c) Average Option Call
tmp = exp(-r*T) * mean(S_aoc, 2) / s0 .* (mean(S_aoc, 2) > K);
delAOC_pathwise = mean(tmp);
gamAOC_pathwise = NaN;

% ----------------------------------
% (4) Likelihood ratio

% (a) Digital Call
s =  tmp1*exp(tmp2*z1);
scoreDC = z1 / (s0*sig*sqrt(T));
del_likelihood = exp(-r*T) * (s > K) .* scoreDC;
delDC_likelihood = mean(del_likelihood);
scoreDC = (z1.^2 - 1 - sig * sqrt(T) * z1) / (s0^2 * sig^2 * T);
gam_likelihood = exp(-r*T) * (s > K) .* scoreDC;
gamDC_likelihood = mean(gam_likelihood);

% (b) Up and Out Call
s =  tmp1*exp(tmp2*z1);
scoreUOC = z1 / (s0*sig*sqrt(T));
del_likelihood = exp(-r*T) * (s < H) .* max(s - K, 0) .* scoreUOC;
delUOC_likelihood = mean(del_likelihood);
scoreUOC = (z1.^2 - 1 - sig * sqrt(T) * z1) / (s0^2 * sig^2 * T);
gam_likelihood = exp(-r*T) * (s < H) .* max(s - K, 0) .* scoreUOC;
gamUOC_likelihood = mean(gam_likelihood);

% (c) Average Option Call
scoreAOC = z_aoc(:, 1) / (s0*sig*sqrt(dt));
del_likelihood = exp(-r*T) * max(mean(S_aoc, 2) - K, 0) .* scoreAOC;
delAOC_likelihood = mean(del_likelihood);
scoreAOC = (z_aoc(:, 1).^2 - 1 - sig * sqrt(dt) * z_aoc(:, 1)) / (s0^2 * sig^2 * dt);
gam_likelihood = exp(-r*T) * max(mean(S_aoc, 2) - K, 0) .* scoreAOC;
gamAOC_likelihood = mean(gam_likelihood);

% ===================================
% ===================================
%
%             Vega
%
% -----------------------------------
tmp1U  =  s0 *exp((r-q-(sig+delSig)^2/2)*T);
tmp1D  =  s0 *exp((r-q-(sig-delSig)^2/2)*T);
tmp2U = (sig+delSig)*sqrt(T);
tmp2D = (sig-delSig)*sqrt(T);

% (1) W/ Common Random Numbers
su = tmp1U*exp(tmp2U*z1);
sd = tmp1D*exp(tmp2D*z1);

% (a) Digital Call
tmpU = exp(-r*T) * (su > K);
tmpD = exp(-r*T) * (sd > K);
vegaDC_CRN = mean((tmpU-tmpD)/(2*delSig));

% (b) Up and Out Call
tmpU = exp(-r*T) * (su < H) .* max(su - K, 0);
tmpD = exp(-r*T) * (sd < H) .* max(sd - K, 0);
vegaUOC_CRN = mean((tmpU-tmpD)/(2*delSig));

% (c) Average Option Call
Su_aoc = zeros(nN,m);
Sd_aoc = zeros(nN,m);
for i = 1:nN
    Su_previous = s0;
    Sd_previous = s0;
    for ite = 1:m
        if ite > 1
            Su_previous = Su_aoc(i, ite - 1);
            Sd_previous = Sd_aoc(i, ite - 1);
        end
        Su_aoc(i,ite) = Su_previous * ...
            exp((r - q - (sig+delSig) * (sig+delSig) / 2) * dt + (sig+delSig) * sqrt(dt) * z_aoc(i, ite));
        Sd_aoc(i,ite) = Sd_previous * ...
            exp((r - q - (sig-delSig) * (sig-delSig) / 2) * dt + (sig-delSig) * sqrt(dt) * z_aoc(i, ite));
    end
end
tmpU = exp(-r*T) * max(mean(Su_aoc, 2) - K, 0);
tmpD = exp(-r*T) * max(mean(Sd_aoc, 2) - K, 0);
vegaAOC_CRN = mean((tmpU-tmpD)/(2*delSig));

% (2) W/O Common Random Numbers

su = tmp1U*exp(tmp2U*z1);
sd = tmp1D*exp(tmp2D*z2);

% (a) Digital Call
tmpU = exp(-r*T) * (su > K);
tmpD = exp(-r*T) * (sd > K);
vegaDC_tilde = mean((tmpU-tmpD)/(2*delSig));

% (b) Up and Out Call
tmpU = exp(-r*T) * (su < H) .* max(su - K, 0);
tmpD = exp(-r*T) * (sd < H) .* max(sd - K, 0);
vegaUOC_tilde = mean((tmpU-tmpD)/(2*delSig));

% (c) Average Option Call
Su_aoc = zeros(nN,m);
Sd_aoc = zeros(nN,m);
for i = 1:nN
    Su_previous = s0;
    Sd_previous = s0;
    for ite = 1:m
        if ite > 1
            Su_previous = Su_aoc(i, ite - 1);
            Sd_previous = Sd_aoc(i, ite - 1);
        end
        Su_aoc(i,ite) = Su_previous * ...
            exp((r - q - (sig+delSig) * (sig+delSig) / 2) * dt + (sig+delSig) * sqrt(dt) * zu_aoc(i, ite));
        Sd_aoc(i,ite) = Sd_previous * ...
            exp((r - q - (sig-delSig) * (sig-delSig) / 2) * dt + (sig-delSig) * sqrt(dt) * zd_aoc(i, ite));
    end
end

tmpU = exp(-r*T) * max(mean(Su_aoc, 2) - K, 0);
tmpD = exp(-r*T) * max(mean(Sd_aoc, 2) - K, 0);
vegaAOC_tilde = mean((tmpU-tmpD)/(2*delSig));

% (3) Pathwise Estimator

% (a) Digital Call
vegaDC_pathwise = NaN;

% (b) Up and Out Call
vegaUOC_pathwise = NaN;

% (c) Average Option Call
tmp = zeros(nN,1);
for i = 1:nN
    factor = 0;
    for j = 1:m
        factor = factor + S_aoc(i, j) * (log(S_aoc(i, j)/s0) - (r - q + sig^2 / 2) * j * dt);
    end 
    tmp(i) = exp(-r*T) * (mean(S_aoc(i, :)) > K) * factor / (m * sig);
end
vegaAOC_pathwise = mean(tmp);

% (4) Likelihood Ratio

% (a) Digital Call
s =  tmp1*exp(tmp2*z1);
scoreDC = 1 / sig * (z1.^2 - sig * sqrt(T) * z1 -1);
vega_likelihood = exp(-r*T) * (s > K) .* scoreDC;
vegaDC_likelihood = mean(vega_likelihood);

% (b) Up and Out Call
s =  tmp1*exp(tmp2*z1);
scoreUOC = 1 / sig * (z1.^2 - sig * sqrt(T) * z1 -1);
vega_likelihood = exp(-r*T) * (s < H) .* max(s - K, 0) .* scoreUOC;
vegaUOC_likelihood = mean(vega_likelihood);

% (c) Average Option Call
scoreAOC = zeros(nN,1);
for i = 1:nN
    tmp = 0;
    for j = 1:m
        tmp = tmp + (z_aoc(i,j)^2 - sig * sqrt(dt) * z_aoc(i,j) - 1);
    end
    scoreAOC(i) = tmp / sig;
end
vega_likelihood = exp(-r*T) * max(mean(S_aoc, 2) - K, 0) .* scoreAOC;
vegaAOC_likelihood = mean(vega_likelihood);

%% Displaying Results

disp('===============');
disp('   Delta');
disp('---------------');
disp('      Digital Call');
disp([' W/ CRN:' num2str(delDC_CRN),  ' W/O CRN:', num2str(delDC_tilde), ...
    ' Pathwise Estimator:', num2str(delDC_pathwise), ' Likelihood Ratio:', num2str(delDC_likelihood)]);
disp(' ');
disp('      Up and Out Call');
disp([' W/ CRN:' num2str(delUOC_CRN),  ' W/O CRN:', num2str(delUOC_tilde), ...
    ' Pathwise Estimator:', num2str(delUOC_pathwise), ' Likelihood Ratio:', num2str(delUOC_likelihood)]);
disp(' ');
disp('      Average Option Call');
disp([' W/ CRN:' num2str(delAOC_CRN),  ' W/O CRN:', num2str(delAOC_tilde), ...
    ' Pathwise Estimator:', num2str(delAOC_pathwise), ' Likelihood Ratio:', num2str(delAOC_likelihood)]);
disp(' ');
disp('===============');
disp('   Gamma');
disp('---------------');
disp('      Digital Call');
disp([' W/ CRN:' num2str(gamDC_CRN),  ' W/O CRN:', num2str(gamDC_tilde), ...
    ' Pathwise Estimator:', num2str(gamDC_pathwise), ' Likelihood Ratio:', num2str(gamDC_likelihood)]);
disp(' ');
disp('      Up and Out Call');
disp([' W/ CRN:' num2str(gamUOC_CRN),  ' W/O CRN:', num2str(gamUOC_tilde), ...
    ' Pathwise Estimator:', num2str(gamUOC_pathwise), ' Likelihood Ratio:', num2str(gamUOC_likelihood)]);
disp(' ');
disp('      Average Option Call');
disp([' W/ CRN:' num2str(gamAOC_CRN),  ' W/O CRN:', num2str(gamAOC_tilde), ...
    ' Pathwise Estimator:', num2str(gamAOC_pathwise), ' Likelihood Ratio:', num2str(gamAOC_likelihood)]);
disp(' ');
disp('===============');
disp('   Vega');
disp('---------------');
disp('---------------');
disp('      Digital Call');
disp([' W/ CRN:' num2str(vegaDC_CRN),  ' W/O CRN:', num2str(vegaDC_tilde), ...
    ' Pathwise Estimator:', num2str(vegaDC_pathwise), ' Likelihood Ratio:', num2str(vegaDC_likelihood)]);
disp(' ');
disp('      Up and Out Call');
disp([' W/ CRN:' num2str(vegaUOC_CRN),  ' W/O CRN:', num2str(vegaUOC_tilde), ...
    ' Pathwise Estimator:', num2str(vegaUOC_pathwise), ' Likelihood Ratio:', num2str(vegaUOC_likelihood)]);
disp(' ');
disp('      Average Option Call');
disp([' W/ CRN:' num2str(vegaAOC_CRN),  ' W/O CRN:', num2str(vegaAOC_tilde), ...
    ' Pathwise Estimator:', num2str(vegaAOC_pathwise), ' Likelihood Ratio:', num2str(vegaAOC_likelihood)]);
disp(' ');