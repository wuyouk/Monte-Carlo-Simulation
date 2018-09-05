% HW3_Q4

% Input

s0 = 200;
T = 1;
sig = 0.35;
r = 0.05;
q = 0.025;
K = 120;

nSims = 20000;
nN = nSims;

% Exact solution
[~, p_ex] = blsprice(s0, K, r, T, sig, q);

% Exact Solution Output
sprintf("Exact price: %.6f", p_ex)

% Standard simulation
z = randn(nSims, 1);
s_st = s0 * exp((r-q-sig*sig/2)*T + sig*sqrt(T)*z);
p_st = exp(-r*T) * max(K - s_st, 0);

% (a) Standard Simulation Output
mean_st = mean(p_st);
sprintf("Standard simulation estimation: %.6f", mean_st)

z = norminv(1 - 0.01 / 2, 0, 1);
st_l = mean_st - z * sqrt(var(p_st)/nN);
st_u = mean_st + z * sqrt(var(p_st)/nN);
sprintf("Standard simulation estimation Confidence Interval: [%.6f, %.6f]", st_l, st_u)

% Stratified Sampling (Sub-optimal)

delta = [-2 -1 1 2];
m = length(delta)+1;
p = zeros(m,1);
phiA = zeros(m,1);
phiB = zeros(m,1);
for i = 1:m
    if i == 1
        p(i) = normcdf(delta(i));
        phiA(i) = 0; phiB(i) = normcdf(delta(i));
    elseif i == m
        p(i) = 1-normcdf(delta(i-1));
        phiA(i) = normcdf(delta(i-1)); phiB(i) = 1;
    else
        p(i) = normcdf(delta(i))-normcdf(delta(i-1));
        phiA(i) = normcdf(delta(i-1)); phiB(i) = normcdf(delta(i));
    end
end

ni = nN * p;

% algorithm
cHat_sub = 0;
sig2N_sub = 0;
for i = 1:m
    n = ceil(ni(i));
    u = rand(n,1);
    u = phiA(i)+(phiB(i)-phiA(i))*u;
    z = norminv(u);
    s_sub = s0 * exp((r-q-sig*sig/2)*T + sig*sqrt(T)*z);
    p_sub = exp(-r*T) * max(K - s_sub, 0);
    
    si = sum(p_sub);
    ui = sum(p_sub.^2);
    theta = si/n;
    sig2 = (ui-si^2/n)/(n-1);
    %
    cHat_sub = cHat_sub + p(i)*theta;
    sig2N_sub = sig2N_sub + sig2*p(i)*p(i)/n;
end

% (b) Sub-optimal Stratified Sampling Estimation Output
mean_sub = cHat_sub;
var_sub = sig2N_sub;
sprintf("Sub-optimal Stratified Sampling: %.6f", mean_sub)

z = norminv(1 - 0.01 / 2, 0, 1);
sub_l = mean_sub - z * sqrt(var_sub);
sub_u = mean_sub + z * sqrt(var_sub);
sprintf("Sub-optimal Stratified Sampling Confidence Interval: [%.6f, %.6f]", sub_l, sub_u)


% Stratified Sampling (Optimal)

% Pilot program to find ni optimal
sig_i = zeros(m,1);
for i=1:m
    np = 5000;
    u = rand(np,1);
    u = phiA(i)+(phiB(i)-phiA(i))*u;
    z = norminv(u);
    s_tmp = s0 * exp((r-q-sig*sig/2)*T + sig*sqrt(T)*z);
    p_tmp = exp(-r*T) * max(K - s_tmp, 0);

    si = sum(p_tmp);
    ui = sum(p_tmp.^2);
    sig_i(i) = sqrt((ui-si^2/np)/(np-1));
end

ns_i = nN*p.*sig_i/sum(p.*sig_i);
% disp([round(ni) round(ns_i)]);

% algorithm
cHat_opt = 0;
sig2N_opt = 0;
for i = 1:m
    n = max(ceil(ns_i(i)),10);
    u=rand(n,1);
    u = phiA(i)+(phiB(i)-phiA(i))*u;
    z = norminv(u);
    
    s_opt = s0 * exp((r-q-sig*sig/2)*T + sig*sqrt(T)*z);
    p_opt = exp(-r*T) * max(K - s_opt, 0);

    si = sum(p_opt);
    ui = sum(p_opt.^2);
    theta = si/n;
    sig2 = (ui-si^2/n)/(n-1);
    %
    cHat_opt = cHat_opt+p(i)*theta;
    sig2N_opt = sig2N_opt+sig2*p(i)*p(i)/n;
end

% (c) Optimal Stratified Sampling Estimation Output
mean_opt = cHat_opt;
var_opt = sig2N_opt;
sprintf("Optimal Stratified Sampling: %.6f", mean_opt)

z = norminv(1 - 0.01 / 2, 0, 1);
opt_l = mean_opt - z * sqrt(var_opt);
opt_u = mean_opt + z * sqrt(var_opt);
sprintf("Optimal Stratified Sampling Confidence Interval: [%.6f, %.6f]", opt_l, opt_u)
