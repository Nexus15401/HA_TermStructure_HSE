function [B] = B_fun(tau,rho,b,mu,sigma,Yt)

mu_tau    = (mu - 0.5*sigma^2)*tau; 
sigma_tau = sigma*sqrt(tau);
    
% A. Distribution parameters (x = epsilon)
    % x ~ N(0,1)
    mu_x    = 0;          % mean vector
    sigma_x = 1;          % variance matrix

% B. Gaussian integration
    nn = 21;                     % order of approximation
    
    % For Equity-TS
    [x,w] = qnwnorm(nn,mu_x,sigma_x);     % Gaussian normal nodes and weights

% C. Function to be integrated
    % For Equity-TS
    f1 = sqrt(1 + b*Yt) - 1;
    f = @(x) f1./sqrt( 1 + b*Yt.*exp(mu_tau + sigma_tau.*x) );

    fexp = w'*f(x);                  % Gaussian integration of f
    fprintf('Guassian Quadrature:      %10.3f\n',fexp)

% F. B (Eq. 57 in P11)
    % B = price of an asset that pays "Y" at the maturity "t+\tau"
    B = (exp(-rho*tau))*fexp;
    