% HW3_Q3

% Input

s0 = 200;
T = 1;
sig = 0.35;
r = 0.05;
q = 0.025;
K = 120;

nSims = 20000;
n = nSims;

% Standard simulation
z = randn(nSims, 1);
s_st = s0 * exp((r-q-sig*sig/2)*T + sig*sqrt(T)*z);
p_st = exp(-r*T) * max(K - s_st, 0);

% Importance Sampling
mu = (log(K / s0) - (r - q - sig * sig /2) / T) / (sig * sqrt(T));
z_is = mu + randn(nSims);
p_is = zeros(nSims,1);
for i = 1:nSims
    s_is = s0 * exp((r-q-sig*sig/2)*T + sig*sqrt(T) .* z_is(i));
    p_is(i) = (exp(-r*T) * max(K - s_is, 0)) * (exp(-mu * z_is(i) + mu ^ 2 / 2));
end
mean(p_is);

% (a) Standard Simulation Output
mean_st = mean(p_st);

sprintf("Standard simulation estimation: %.6f", mean_st)

z = norminv(1 - 0.01 / 2, 0, 1);
st_l = mean_st - z * sqrt(var(p_st)/n);
st_u = mean_st + z * sqrt(var(p_st)/n);
sprintf("Standard simulation estimation Confidence Interval: [%.6f, %.6f]", st_l, st_u)

% (b) Importance Sampling Output
mean_is = mean(p_is);

sprintf("Importance Sampling Estimation: %.6f", mean_is)

z = norminv(1 - 0.01 / 2, 0, 1);
is_l = mean_is - z * sqrt(var(p_is)/n);
is_u = mean_is + z * sqrt(var(p_is)/n);
sprintf("Importance Sampling Estimation Confidence Interval: [%.6f, %.6f]", is_l, is_u)
