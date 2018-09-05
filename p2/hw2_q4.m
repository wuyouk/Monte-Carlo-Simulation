% Input
r = 0.04;
q = 0.015;
sig = 0.30;
T = 1;
maturity = T;
Spot = 1000;
K = 1100;

% Simulation
m = 12;
dt = maturity/m;

% Pilot to find c
n = 30;
S = zeros(n,m);
Z1 = zeros(n,1);
Z2 = zeros(n,1);
Z3 = zeros(n,1);

for i = 1:n
    S(i, 1) = Spot;
    S_previous = Spot;
    for ite = 1:m
        if ite > 1
            S_previous = S(i, ite - 1);
        end
        z = randn;
        S(i, ite) = S_previous * exp((r - q - sig * sig / 2) * dt + sig * sqrt(dt) * z);
    end
    Z1(i) = S(i, m);
    Z2(i) = exp(-r * T) * max(S(i, m) - K, 0);
    Z3(i) = exp(-r * T) * max(geomean(S(i, :)) - K, 0);
end
Y = exp(-r * T) * max(mean(S,2) - K, 0);
E1 = Spot * exp((r-q) * T);
[E2, P2] = blsprice(Spot, K, r, maturity, sig, q);
sig_new = sig / m * sqrt((m+1)*(2*m+1)/6);
q_new = -((m+1)/(2*m)*(r - q - sig^2 / 2) + 1/2*sig_new^2 - r);
[E3, P3] = blsprice(Spot, K, r, maturity, sig_new, q_new);
c1 = pilot_c(Y, Z1, E1, n);
c2 = pilot_c(Y, Z2, E2, n);
c3 = pilot_c(Y, Z3, E3, n);

% Main simulation with Control Variates
n = 50000;
S = zeros(n,m);
Z1 = zeros(n,1);
V1 = zeros(n,1);
Z2 = zeros(n,1);
V2 = zeros(n,1);
Z3 = zeros(n,1);
V3 = zeros(n,1);

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
    Yi = exp(-r * T) * max(mean(S(i, :)) - K, 0);
    Z1(i) = S(i, m);
    V1(i) = Yi + c1 * (Z1(i) - E1);
    Z2(i) = exp(-r * T) * max(S(i, m) - K, 0);
    V2(i) = Yi + c2 * (Z2(i) - E2);
    Z3(i) = exp(-r * T) * max(geomean(S(i, :)) - K, 0);
    V3(i) = Yi + c3 * (Z3(i) - E3);
end

% Standard Result
prices_standard = max(mean(S,2) - K, 0) * exp(-r * T);
sol_standard = mean(prices_standard);
sprintf("Standard simulation price is %.3f", sol_standard)

% Standard Confidence Interval
z = norminv(1 - 0.01 / 2, 0, 1);
cil_CLT = sol_standard - z * sqrt(var(prices_standard)/n);
ciu_CLT = sol_standard + z * sqrt(var(prices_standard)/n);
sprintf("Standard simulation Confidence Interval is [%.3f, %.3f]", cil_CLT, ciu_CLT)

% Output results with Control Variates
output('b', mean(V1), var(V1), n);
output('c', mean(V2), var(V2), n);
output('d', mean(V3), var(V3), n);

% Output function for theta and confidence interval
function output(index, theta, variance, n)
sprintf("(%s) Simulation with Control Variates price is %.3f", index, theta)
z = norminv(1 - 0.01 / 2, 0, 1);
cil_cvCLT = theta - z * sqrt(variance/n);
ciu_cvCLT = theta + z * sqrt(variance/n);
sprintf("(%s) Simulation with Control Variates Confidence Interval is [%.3f, %.3f]", ...
    index, cil_cvCLT, ciu_cvCLT)
end

% Pilot program to find C*
function [c] = pilot_c(ys, zs, ez, n)
% covariance = cov(ys, zs);
% covYZ = covariance(1,2);
covYZ = sum((ys-mean(ys)).*(zs-ez))/(n-1);
varZ = sum((zs - ez).^2)/(n-1);
c = -covYZ/varZ;
end