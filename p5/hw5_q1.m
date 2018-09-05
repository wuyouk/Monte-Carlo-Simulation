function hw5_q1(mu, sig, n1, caseNumber)

%% ============================================================
%   Conjugate Prior for the Normal Distribution
%    
%   Find Conjugate prior distribution for \mu & \sigma^2
%  -------------------------------------------------------------


%% Likelihood

%mu  = 4.0;
%sig = 2.5;
%n1 = 100;

xRange = -10:.002:10;
yRange = normpdf(xRange, mu, sig);

plot(xRange, yRange, 'k'); 

y = normrnd(mu, sig, n1, 1);
yHat = mean(y);




%%
if strcmp(caseNumber,'1') == 1
    
    %% Case 1: \sigma^2 is known, \mu is random
    
    mu0 = -3;
    sig0 = 1;
    yPrior = normpdf(xRange, mu0, sig0);
    plot(xRange, yPrior, 'r');
    
    for j = 1:n1
        
        %posterior
        eta = (sig^2*mu0 + y(j)*sig0^2)/(sig^2 + sig0^2);
        kappa = sqrt((sig^2 * sig0^2)/(sig^2 + sig0^2));
        
        %disp([eta kappa]);
        
        yPosterior = normpdf(xRange, eta, kappa);
        plot(xRange, yPosterior, 'g');
        hold on;
        yMax = max(yPosterior);
        plot(mu, 0, 'o', 'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
        plot([mu, mu], [0,yMax], 'r');
        title('Case I -- \mu_{posterior} (\sigma^2 is known, \mu is unknown)');
        hold off;
        
        pause;
        
        %update prior by using posterior as a prior
        mu0 = eta;
        sig0 = kappa;
        
    end
    
elseif strcmp(caseNumber,'2') == 1
        
    %% Case 2: \mu is known, \sigma^2 is random
    a0 = 0.2;
    b0 = 0.5;
    lRange = 0.002:.002:0.5;
    
    yPrior = gampdf(lRange, a0, 1/b0);
    plot(lRange, yPrior, 'r');
    
    
    for j = 1:n1
        
        %posterior
        aNew = a0 + 1/2;
        bNew = b0 + sum((y(j)-mu).^2)/2;
        
        %disp([aNew bNew]);
        
        yPosterior = gampdf(lRange, aNew, 1/bNew);
        plot(lRange, yPosterior, 'g');
        
        %update prior by using posterior parameters for the prior
        a0 = aNew;
        b0 = bNew;
        
        plot(lRange, yPosterior);
        hold on;
        yMax = max(yPosterior);
        plot(1/sig^2,0, 'o', 'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
        plot([1/sig^2, 1/sig^2], [0,yMax], 'r');
        title('Case II -- \sigma^2_{posterior} (\sigma^2 is unknown, \mu is known)');
        hold off;
        
        pause;
        
    end
    

elseif strcmp(caseNumber,'3') == 1
    
    %% Case 3: mu and sig^2 are unknown
    % To be extended
    
    ldaRange = 0.002:.002:3;
    muRange = 0.01:0.01:10;
    
    a0 = 0.2;
    b0 = 0.5;
    
    k0 = 1/(1.0^2);
    mu0 = 0.5;
    
    ldaPrior = gampdf(ldaRange, a0, 1/b0);
    
    t0 = (muRange-mu0)/sqrt(b0/(a0*k0));
 
    muPrior = tpdf(t0, 2*a0);
    subplot(2,1,1);
    plot(ldaRange,ldaPrior,'r');
    subplot(2,1,2);
    plot(muRange,muPrior,'r');
    
    pause;
    
    
    
    for j = 1:n1
        
        % lambda
        muN = (k0*mu0+1*y(j))/(k0+1);
        kN = k0+1;

        aN = a0 + 1/2;
        bN = b0 + k0*1*(y(j)-mu0)^2/(2*(k0+1));

        ldaPosterior = gampdf(ldaRange, aN, 1/bN);
        subplot(2,1,1);
        plot(ldaRange, ldaPosterior, 'g');
        hold on;
        yMax = max(ldaPosterior);
        plot(1/sig^2,0, 'o', 'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
        plot([1/sig^2, 1/sig^2], [0,yMax], 'r');
        title('Case III -- \sigma^2_{posterior} (Both \sigma^2 & \mu are unknown) ');
        hold off;
        
        % mu
        t = (muRange-muN)/sqrt(bN/(aN*kN));
        muPosterior = tpdf(t, 2*aN);

        subplot(2,1,2);
        plot(muRange, muPosterior, 'g');
        hold on;
        yMax = max(muPosterior);
        plot(mu,0, 'o', 'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
        plot([mu, mu], [0,yMax], 'r');
        title('Case III -- \mu_{posterior} (Both \sigma^2 & \mu are unknown)');
        hold off;
        
        % update
        a0 = aN;
        b0 = bN;
        mu0 = muN;
        k0 = kN;
        
        pause;
        
        
    end
    
    
    
    
    
    

    
end
