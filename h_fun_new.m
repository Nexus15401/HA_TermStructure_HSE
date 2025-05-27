function [h,rtau] = h_fun_new(tau,rho,lambda,mu,sigma,gamma1,gamma2,...
                         lambdax,xt,wt,tildec1t)
%{
P19
Calculate Price (h) and Equity Yields (rtau)
Hamilton Galindo
May 2024
%}

% Price ==> h = A*E[f]

% A. Preliminary
A    = ( exp(-rho*tau) )*( exp(xt + gamma1*wt) )*tildec1t^(gamma1);
wbar = (mu - (1/2)*sigma^2)/lambdax;
a = (1-lambda)/lambda;

mu_tau    = wbar + (wt - wbar)*exp(-lambdax*tau); 
sigma_tau = ( (sigma^2)/(2*lambdax) )*( 1 - exp(-2*lambdax*tau));
    
% B. Distribution parameters
    % w_{t+tau} = x ~ N(mu_tau,sigma^2_tau)
    mu_x    = mu_tau;          % mean vector
    sigma_x = sigma_tau;          % variance matrix

% C. Gaussian integration
    nn = 21;                     % order of approximation
    
    % For Equity-TS
    [xx,weights] = qnwnorm(nn,mu_x,sigma_x);     % Gaussian normal nodes and weights
                                                 % xx = w_{t+tau}  
    % Consumption at "t+tau" (for every "xx" value)
    c10 = 0.5; %initial point of agent1's consumption
    for i=1:nn
        c1_ttau(i) = fsolve(@(c1) c_op(c1,a,gamma1,gamma2,xx(i)),c10);
    end
    
    % Function to be integrated
    f = ( exp( -(gamma1 - 1)*xx ) ).*( c1_ttau'.^(-gamma1) );

    Ef = weights'*f;                  % Gaussian integration of f
    fprintf('Guassian Quadrature:      %10.3f\n',Ef)

% D. Price = h
    % h = price of an asset that pays "Y" at the maturity "t+\tau"
    h = A*Ef;
    
% E. Equity yield
B = ( exp( rho*tau - (gamma1-1)*wt) )*tildec1t^(-gamma1);
rtau = (1/tau)*log(B/Ef);