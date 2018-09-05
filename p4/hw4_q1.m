% Pricing Heston Utilizing various schemes

% Input

s0 = 100;
K = 110;
T = 1;
sig = 0.30;
r = 0.025;
q = 0.0125;
theta = 0.0625;
ki = 2.75;
lambda = 0.0125;
v0 = 0.05;
rho = -0.65;

[c, p] = blsprice(s0, K, r, T, sig, q);

payOff1 = 0;
payOff2 = 0;
payOff3 = 0;

nSims = 1000000;
m  = 52;
dt = T/m;
for j = 1:nSims
    s1 = s0;
    s2 = s0;
    s3 = s0;
    v1 = v0;
    v2 = v0;
    v3 = v0;
    for i = 1:m
        z1 = randn;
        z3 = randn;
        z2 = rho * z1 + sqrt(1 - rho^2) * z3;
        % Euler
        s1 = s1 + (r - q) * s1 * dt + sqrt(v1) * s1 * sqrt(dt) * z1;
        v1 = v1 + ki * (theta - v1) * dt + lambda * sqrt(v1) * sqrt(dt) * z2;
        % Milstein
        s2 = s2 + (r - q) * s2 * dt + sqrt(v2) * s2 * sqrt(dt) * z1 + 0.5 * v2 * s2 * dt * (z1^2 - 1);
        v2 = v2 + ki * (theta - v2) * dt + lambda * sqrt(v2) * sqrt(dt) * z2 + 0.5^2 * lambda^2 * dt * (z2^2 - 1);
        % Runge-Kutta
        s3_tilde = s3 + (r - q) * s3 * dt + sqrt(v3) * s3 * sqrt(dt);
        s3 = s3 + (r - q) * s3 * dt + sqrt(v3) * s3 * sqrt(dt) * z1 + ...
            1 / (2 * sqrt(dt)) * sqrt(v3) * (s3_tilde - s3) * dt * (z1^2 - 1);
        v3_tilde = v3 + ki * (theta - v3) * dt + lambda * sqrt(v3) * sqrt(dt);
        v3 = v3 + ki * (theta - v3) * dt + lambda * sqrt(v3) * sqrt(dt) * z2 + ...
            1 / (2 * sqrt(dt)) * lambda * (sqrt(v3_tilde) - sqrt(v3)) * dt * (z2^2 - 1);   
    end
    payOff1 = payOff1 + max(s1-K,0);
    payOff2 = payOff2 + max(s2-K,0);
    payOff3 = payOff3 + max(s3-K,0);

end

c1 = exp(-r*T)*payOff1/nSims;
c2 = exp(-r*T)*payOff2/nSims;
c3 = exp(-r*T)*payOff3/nSims;



sprintf(" Euler: %.6f \n Milstein: %.6f \n Runge-Kutta: %.6f \n", c1, c2, c3)

%% Parameters

alpha = 1.5;
eta = 0.1;
n = 10;
N = 2^n;

S0 = 100;
K = 110;
r = 0.025;
q = 0.0125;
T = 1.0;

%% Model

model = 'Heston';

switch model
    
    case 'GBM'
        
        sig = 0.3;
        params = sig;
        
        call_GBM = blsprice(S0, K, r, T, sig, q);
        call = genericFFT(model, eta, alpha, N, S0, K, r, q, T, params);
        
        disp([call_GBM call]);
        
    case 'Heston'
        
        kappa = 2.75;
        theta= 0.0625;
        lambda = 0.0125;
        rho = -0.65;
        v_0 = 0.05;
        
        params = [kappa theta lambda rho v_0];
        call = genericFFT(model, eta, alpha, N, S0, K, r, q, T, params);
        
end

sprintf(" Exact solution: %.6f \n", call)

function call = genericFFT(model, eta, alpha, N, S0, K, r, q, T, params)

lda = 2*pi/(N*eta);
%
j = 0:N-1;
nu = j*eta; 
nu = nu(:);

C = exp(-r*T); % constant

CF_u = modelCF(model, nu-(alpha+1)*1i, S0, r, q, T, params);

beta = log(K);

FourierTransformOfModifiedCall = zeros(N,1);	
for i = 1:N
    if i == 1
        FourierTransformOfModifiedCall(i) = (eta/2) * C * exp(-1i*beta*nu(i)) * CF_u(i) / ( (alpha + 1i*nu(i)) * (alpha + 1i*nu(i) + 1) );
    else
        FourierTransformOfModifiedCall(i) = eta     * C * exp(-1i*beta*nu(i)) * CF_u(i) / ( (alpha + 1i*nu(i)) * (alpha + 1i*nu(i) + 1) );
    end
end

y = fft(FourierTransformOfModifiedCall);

a = zeros(N,1);
for m = 1:N
	a(m) = exp(-alpha*(beta+(m-1)*lda))/pi;
end

c = (real(y)).*a;

call = c(1);
end

function cf = modelCF(model, u, S0, r, q, T, params)


switch model
    
    case 'GBM'
        
        sig = params(1);
        cf = exp( 1i*(log(S0)+(r-q-sig*sig/2)*T)*u - 0.5*sig*sig*u.*u*T );
        
    case {'Heston', 'GBMSA'}
        
        %params = [kappa, theta, sigma, rho, v_0];
        kappa = params(1);
        theta = params(2);
        sigma = params(3);
        rho   = params(4);
        v_0   = params(5);
        
        g = sqrt((kappa-1i*rho*sigma.*u).^2+(u.*u+1i.*u)*sigma*sigma);
        beta = kappa-rho*sigma*1i.*u;
        tmp = g*T/2;
        %
        temp1 = 1i*(log(S0)+(r-q)*T)*u + kappa*theta*T*beta/(sigma*sigma);
        temp2 = -(u.*u+1i*u)*v_0./(g.*coth(tmp)+beta);
        temp3 = (2*kappa*theta/(sigma*sigma)).*log(cosh(tmp)+(beta./g).*sinh(tmp));
        
        cf = exp(temp1+temp2-temp3);
        
        
    case 'VG'
        
    case 'Merton'
        
    case 'VGSA'
        
    case 'VGSSD'
        
    case 'CGMY'
        
    case 'CGMYSA'
        
end
end