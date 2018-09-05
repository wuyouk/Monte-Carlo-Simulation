function hw5_q3(theta, v, a0, b0, n1)

%% ============================================================
%   Conjugate Prior for the Gamma Distribution
%  -------------------------------------------------------------


%% Likelihood
%a0  = 3.0;
%b0 = 3.0;
%n1 = 100;

xRange = 0.002:.002:5;
yRange = gampdf(xRange, a0, 1/b0);

plot(xRange, yRange, 'k'); 

y = gamrnd(v, 1/theta, n1, 1);

for j = 1:n1
        
    %posterior
	aNew = a0 + 1*v;
	bNew = y(j) + b0;
        
    %disp([aNew bNew]);
        
    yPosterior = gampdf(xRange, aNew, 1/bNew);
	plot(xRange, yPosterior, 'g');
        
	%update prior by using posterior parameters for the prior
	a0 = aNew;
	b0 = bNew;
        
	hold on;
	yMax = max(yPosterior);
	plot(theta,0, 'o', 'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
	plot([theta, theta], [0,yMax], 'r');
	title('\theta_{posterior}');
	hold off;
        
	pause;
        
end

