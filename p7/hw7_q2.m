%% Parameters
r = 0.10; q = 0.01; s0 = 100; K = 100; sig = 0.40; tau = 1;

m = 365;
numPaths = 500000;

dt = tau/m;
indicators = zeros(numPaths,1);

%% Initiate paths (not efficient, very expensive)
% should utilize the Brownian Bridge
%{
s = zeros(m+1,numPaths);
t = zeros(m+1,numPaths);
for j = 1:numPaths
    s(1,j) = s0;
    t(1,j) = 0;
    T = tau;
    for i = 2:m+1
        T = T-dt;
        z = randn;
        s(i,j) = s(i-1,j) * exp((r-q-sig*sig/2)*dt + sig*sqrt(dt)*z);
        t(i,j) = tau - T;
    end
end
%}
sMin = 10;
sMax = 400;
% 
typeOfPolynomial = 'Chebychev_secondKind';

%% Algorithm
warning off;

%P = max(K-s(m+1,:),0);


vHat = zeros(numPaths,1);

s = zeros(m+1,numPaths);
t = zeros(m+1,numPaths);

T = tau;

for i = m:-1:2
    
    if i == m
        z = randn(numPaths,1);
        WT = sqrt(T)*z;
        s(m+1,:) = s0 * exp((r-q-sig*sig/2)*T + sig*sqrt(T)*z);
        t(m+1,:) = T;
        P = max(K-s(m+1,:),0);
    end

    z = randn(numPaths,1);
	time = T - dt;
	Wt = (time / T) .* WT + sqrt((T - time) * time / T).*z;
	s(i,:) = s0 * exp((r-q-sig*sig/2)*time + sig*Wt);
	t(i,:) = time;
    
    % update T and WT to current point
	T = time;
	WT = Wt;
    
    s_i = s(i,:);
    g = max(K-s_i,0);
    
    % in-the-money indicator
    indicator = (g>0);
    disp(sum(indicator));
    
    xi = s_i(indicator);
    xi = xi(:);
    
    if strcmp(typeOfPolynomial,'Laguerre') == 1
	
        xXi = hw7_q2_constructX(xi, typeOfPolynomial);
		
    elseif (strcmp(typeOfPolynomial,'Chebychev_firstKind') == 1) || (strcmp(typeOfPolynomial,'Chebychev_secondKind') == 1)
	
        % make it shifted from (-1,1) to (sMin, sMax)
        xXi = hw7_q2_constructX(2*(xi-sMin)/(sMax-sMin)-1, typeOfPolynomial);
        
    end
    
    yi = exp(-r*dt)*P(indicator);
    yi = yi(:);

    % regression
    alpha = regress(yi, xXi);
	
    vH = xXi*alpha;
    vHat(indicator==1) = vH;
    
    for j = 1:numPaths
        if indicator(j) == 1 && g(j) > vHat(j)
            P(j) = g(j);
        else
            P(j) = exp(-r*dt)*P(j);
        end
    end
end

premiumHat = exp(-r*dt)*mean(P);
disp(premiumHat);

%%