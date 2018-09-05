% HW3_Q2

% Input

s0 = 100;
T = 1;
t1 = T / 2;
sig = 0.30;
r = 0.04;
q = 0.015;
H = 110;
l1 = 0.9;
l2 = 1.1;

nSims = 20000;
n = nSims;
discountedPayoffs_st = zeros(nSims,1);
discountedPayoffs_cmc = zeros(nSims,1);

for i = 1:nSims
    z1 = randn;
    s1 = s0 * exp((r-q-sig*sig/2)*t1 + sig*sqrt(t1)*z1);

    % Usual simulation
    z2 = randn;
    s2 = s1 * exp((r-q-sig*sig/2)*(T-t1) + sig*sqrt(T-t1)*z2);
    if s1 < H
        discountedPayoffs_st(i) = exp(-r*T) * max(l1 * s1 - s2,0);
    else
        discountedPayoffs_st(i) = exp(-r*T) * max(s2 - l2 * s1,0);
    end
    
    % Conditional MC
    [~, p_cmc] = blsprice(s1, l1 * s1, r, T-t1, sig, q);
    [c_cmc, ~] = blsprice(s1, l2 * s1, r, T-t1, sig, q);
    if s1 < H
        discountedPayoffs_cmc(i) = exp(-r*t1) * p_cmc;
    else
        discountedPayoffs_cmc(i) = exp(-r*t1) * c_cmc;
    end
end

% (a) Standard Simulation Output
mean_st = mean(discountedPayoffs_st);

sprintf("Standard simulation estimation: %.6f", mean_st)

z = norminv(1 - 0.01 / 2, 0, 1);
st_l = mean_st - z * sqrt(var(discountedPayoffs_st)/n);
st_u = mean_st + z * sqrt(var(discountedPayoffs_st)/n);
sprintf("Standard simulation estimation Confidence Interval: [%.6f, %.6f]", st_l, st_u)

% (b) Conditional MC Output
mean_v = mean(discountedPayoffs_cmc);

sprintf("Conditional Monte Carlo simulation estimation: %.6f", mean_v)

z = norminv(1 - 0.01 / 2, 0, 1);
v_l = mean_v - z * sqrt(var(discountedPayoffs_cmc)/n);
v_u = mean_v + z * sqrt(var(discountedPayoffs_cmc)/n);
sprintf("Conditional Monte Carlo simulation estimation Confidence Interval: [%.6f, %.6f]", v_l, v_u)
