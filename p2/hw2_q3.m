% Input
r = 0.04;
q = 0.015;
sig = 0.30;
T = 1;
maturity = T;
Spot = 1000;
K = 1100;

% Standard simulation
m = 12;
n = 50000;
dt = maturity/m;

S = zeros(n,m);

for i = 1:n
    S(i, 1) = Spot;
    S_previous = Spot;
    for ite = 1:12
        if ite > 1
            S_previous = S(i, ite - 1);
        end
        z = randn;
        S(i, ite) = S_previous * exp((r - q - sig * sig / 2) * dt + sig * sqrt(dt) * z);
    end
end

prices_standard = max(mean(S,2) - K, 0) * exp(-r * T);
sol_standard = mean(prices_standard);
sprintf("Standard simulation price is %.3f", sol_standard)

% Standard Confidence Interval
z = norminv(1 - 0.01 / 2, 0, 1);
cil_CLT = sol_standard - z * sqrt(var(prices_standard)/n);
ciu_CLT = sol_standard + z * sqrt(var(prices_standard)/n);
sprintf("Standard simulation Confidence Interval is [%.3f, %.3f]", cil_CLT, ciu_CLT)

% Antithetic Variates
n = 50000/2;
dt = maturity/m;
S1 = zeros(n,m);
S2 = zeros(n,m);
for i = 1:n
    S1(i, 1) = Spot;
    S2(i, 1) = Spot;
    S1_previous = Spot;
    S2_previous = Spot;
    for ite = 1:12
        if ite > 1
            S1_previous = S1(i, ite - 1);
            S2_previous = S2(i, ite - 1);
        end
        z = randn;
        S1(i, ite) = S1_previous * exp((r - q - sig * sig / 2) * dt + sig * sqrt(dt) * z);
        S2(i, ite) = S2_previous * exp((r - q - sig * sig / 2) * dt - sig * sqrt(dt) * z);
    end
end

prices_av = (max(mean(S1,2) - K, 0) + max(mean(S2,2) - K, 0))/2 * exp(-r * T);
sol_av = mean(prices_av);
sprintf("Simulation with antithetic variates price is %.3f", sol_av)

% Simulation with antithetic variates Confidence Interval
z = norminv(1 - 0.01 / 2, 0, 1);
cil_avCLT = sol_av - z * sqrt(var(prices_av)/n);
ciu_avCLT = sol_av + z * sqrt(var(prices_av)/n);
sprintf("Simulation with antithetic variates Confidence Interval is [%.3f, %.3f]", cil_avCLT, ciu_avCLT)
