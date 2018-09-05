% Input
r = 0.05;
q = 0.01;
sig = 0.30;
T = 1;
maturity = T;
Spot = 100;
K = 90;

% Standard simulation
m = 250;
n = 10000;
dt = maturity/m;
S = zeros(n,1);

for i = 1:n
    S(i) = Spot;
    for ite = 1:m
        z = randn;
        S(i) = S(i) * exp((r - q - sig * sig / 2) * dt + sig * sqrt(dt) * z);
    end
end

prices_standard = max(S - K, 0) * exp(-r * T);
sol_standard = mean(prices_standard);
sprintf("Standard simulation price is %.3f", sol_standard)

% Exact solution
[C, P] = blsprice(Spot, K, r, maturity, sig, q);
sol_exact = C;
sprintf("Exact solution price is %.3f", sol_exact)

% CLT Confidence Interval
z = norminv(1 - 0.025, 0, 1);
cil_CLT = sol_standard - z * sqrt(var(prices_standard)/n);
ciu_CLT = sol_standard + z * sqrt(var(prices_standard)/n);
sprintf("Confidence Interval based on Central Limit Theorem is [%.3f, %.3f]", cil_CLT, ciu_CLT)

% Bootstrap Con?dence Interval
% # of bootstrap samples
B = 100;
V_B = zeros(B,1);
[f1, x1] = myEmpiricalCDF(prices_standard);

thetaHat = mean(prices_standard);
for l = 1:B 
    uU = rand(n,1);
    xX = interp1(f1, x1, uU, 'linear', 'extrap');
    theta_B = mean(xX);
	V_B(l) = (theta_B - thetaHat)^2;
end
mse_F = mean(V_B);
z = norminv(1 - 0.05 / 2, 0, 1);

cil_BSP = thetaHat - z * sqrt(mse_F);
ciu_BSP = thetaHat + z * sqrt(mse_F);
sprintf("Confidence Interval based on Bootstrap is [%.3f, %.3f]", cil_BSP, ciu_BSP)

% Construct Empirical CDF
function [f, x] = myEmpiricalCDF(x)
tinyNumber = 0.00000001;
x_sorted = sort(x, 'ascend');
len2 = length(x);
freq = ones(len2,1)/len2;
f = [0; cumsum(freq)];
x = [x_sorted(1)-tinyNumber; x_sorted];
end
