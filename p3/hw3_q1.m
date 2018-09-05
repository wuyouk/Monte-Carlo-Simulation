% HW3_Q1
% (a) Standard simulation
n = 10000;

u1 = rand(n,1);
u2 = rand(n,1);
l1 = 2;
l2 = 1;
x1 = -1 / l1 * log(1 - u1);
x2 = -1 / l2 * log(1 - u2);
x = (x1 + x2 > 8);
mean_x = mean(x);
sprintf("Standard simulation estimation: %.6f", mean_x)

z = norminv(1 - 0.01 / 2, 0, 1);
x_l = mean_x - z * sqrt(var(x)/n);
x_u = mean_x + z * sqrt(var(x)/n);
sprintf("Standard simulation estimation Confidence Interval: [%.6f, %.6f]", x_l, x_u)

% (b) Conditional MC
v = (x1 <= 8) .* exp(x1 - 8) + (x1 > 8) .* 1;
mean_v = mean(v);
sprintf("Conditional Monte Carlo simulation estimation: %.6f", mean_v)

z = norminv(1 - 0.01 / 2, 0, 1);
v_l = mean_v - z * sqrt(var(v)/n);
v_u = mean_v + z * sqrt(var(v)/n);
sprintf("Conditional Monte Carlo simulation estimation Confidence Interval: [%.6f, %.6f]", v_l, v_u)
