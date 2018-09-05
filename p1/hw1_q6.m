Spot = 100;
K = 90;
rfr = 0.05;
maturity = 1;
sig = 0.35;
q = 0.07;
[C, P] = blsprice(Spot, K, rfr, maturity, sig, q);

m = 52;
dt = maturity/m;
nSim = 1000;

syntheticC = zeros(nSim,1);
for j = 1:nSim
    S = Spot;
    T = maturity;
    SPrevious = 0;
    delCPrevious = 0;
    for i = 1:m
        [delC, delP] = blsdelta(S, K, rfr, T, sig, q);
        
        %disp([S T delC]);
        %pause(1);
        
        syntheticC(j) = syntheticC(j) + exp(rfr*(m-i)*dt) * ((delC-delCPrevious)*S - (SPrevious * delCPrevious * (exp(q * dt) - 1)));
        T = T-dt;
        z = randn;
        SPrevious = S;
        S = SPrevious * exp((rfr-q-sig*sig/2)*dt + sig*sqrt(dt)*z);
        delCPrevious = delC;
    end
    syntheticC(j) = syntheticC(j) - delC*S + max(S-K,0);
    
end
%plot(syntheticC);
disp([C mean(exp(-rfr*maturity)*syntheticC) std(exp(-rfr*maturity)*syntheticC)]);