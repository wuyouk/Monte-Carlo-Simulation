u = rand(1000,1);

theta = mean(1/sqrt(2*pi)*exp(-(u.^2)/2)) + 0.5;

disp(theta)

