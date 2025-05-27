function stochastic_bond_pricing()
    % Parameters
    rho = 0.05;         % Discount rate
    mu = 0.01;          % Drift coefficient
    sigma = 0.1;        % Volatility
    Y0 = 1;             % Initial value
    T = 10;             % Total time
    dt = 0.01;          % Time step
    N = T/dt;           % Number of steps
    M = 1000;           % Number of simulations
    tau = 1;            % Time horizon for bond pricing
    lambda = 0.5;       % Investor weight
    gamma1 = 1;         % Risk aversion parameter (assuming CRRA utility)
    
    % Derived parameters
    b = 4 * ((1-lambda)/lambda)^(2/gamma1);
    
    % Time vector
    t = 0:dt:T;
    
    % Initialize arrays
    Y = zeros(M, length(t));
    Y(:,1) = Y0;
    
    % Simulate M paths of Y_t
    for i = 1:M
        dW = sqrt(dt) * randn(1, N);
        for j = 2:length(t)
            Y(i,j) = Y(i,j-1) * exp((mu - 0.5*sigma^2)*dt + sigma*dW(j-1));
        end
    end
    
    % Calculate r_t = (sqrt(1+b*Y_t)-1)
    r = sqrt(1 + b*Y) - 1;
    
    % Calculate B_t = E[r_t / (e^{-rho*tau} * r_{t+tau})]
    B = zeros(1, length(t)-tau/dt);
    for k = 1:length(t)-tau/dt
        current_idx = k;
        future_idx = k + tau/dt;
        B(k) = mean(r(:,current_idx) ./ (exp(-rho*tau) * r(:,future_idx)));
    end
    
    % Calculate consumption paths (c1 and c2)
    c1 = (2/b) * (sqrt(1 + b*Y) - 1);
    c2 = Y - c1;
    
    % Visualization
    figure('Position', [100, 100, 1200, 800]);
    
    % Plot sample paths of Y_t
    subplot(2,2,1);
    plot(t, Y(1:5,:));
    title('Sample Paths of Y_t');
    xlabel('Time');
    ylabel('Y_t');
    grid on;
    
    % Plot mean and quantiles of Y_t
    subplot(2,2,2);
    meanY = mean(Y);
    quantY = quantile(Y, [0.05, 0.95]);
    plot(t, meanY, 'b', 'LineWidth', 2);
    hold on;
    plot(t, quantY(1,:), 'r--');
    plot(t, quantY(2,:), 'r--');
    title('Mean and 90% Confidence Bands of Y_t');
    xlabel('Time');
    ylabel('Y_t');
    legend('Mean', '5% Quantile', '95% Quantile');
    grid on;
    
    % Plot bond price B_t
    subplot(2,2,3);
    plot(t(1:end-tau/dt), B);
    title('Bond Price B_t');
    xlabel('Time');
    ylabel('B_t');
    grid on;
    
    % Plot consumption paths
    subplot(2,2,4);
    plot(t, mean(c1), 'b', t, mean(c2), 'r');
    title('Average Consumption Paths');
    xlabel('Time');
    ylabel('Consumption');
    legend('Investor 1 (c_{1t})', 'Investor 2 (c_{2t})');
    grid on;
    
    % Additional analysis: Histogram of Y_T
    figure;
    histogram(Y(:,end), 50);
    title('Distribution of Y_T at Final Time');
    xlabel('Y_T');
    ylabel('Frequency');
    
    % Save workspace variables
    save('stochastic_results.mat', 't', 'Y', 'B', 'c1', 'c2');
end