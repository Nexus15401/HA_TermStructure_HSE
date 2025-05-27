function stochastic_bond_pricing_adapted()
    %% P1 - Adapted Bond Pricing Version
    % Author: [Your Name]
    % Based on: Hamilton Galindo's two-agent model framework
    % Adapted for bond pricing with stochastic processes
    %==========================================================================
    clear; close all;

    %% Model Calibrations (Three Cases)
    names = {'Case 1: Baseline', 'Case 2: High Volatility', 'Case 3: Low Discount'};
    lambda = [0.5, 0.6, 0.4];      % Weight parameters
    rho = [0.05, 0.06, 0.03];      % Discount rates
    mu = [0.01, 0.015, 0.008];     % Drift coefficients
    sigma = [0.1, 0.2, 0.05];      % Volatility parameters
    gamma1 = [2, 5, 1];            % Risk aversion coefficients

    %% Simulation Parameters
    Y0 = 1;             % Initial endowment
    T = 10;             % Time horizon
    dt = 0.01;          % Time step
    N = T/dt;           % Number of steps
    M = 1000;           % Number of simulations
    tau = 1;            % Bond maturity

    %% Preallocate Results Structure
    results = cell(1,3);

    %% Main Simulation Loop for Each Calibration
    for i = 1:3
        % Derived parameters
        b = 4 * ((1-lambda(i))/lambda(i)^(2/gamma1(i));
        
        % Initialize arrays
        Y = zeros(M, N+1);
        Y(:,1) = Y0;
        t = 0:dt:T;
        
        % Simulate M paths of Y_t (Geometric Brownian Motion)
        for m = 1:M
            dW = sqrt(dt) * randn(1, N);
            for j = 2:N+1
                Y(m,j) = Y(m,j-1) * exp((mu(i) - 0.5*sigma(i)^2)*dt + sigma(i)*dW(j-1));
            end
        end
        
        % Calculate bond pricing components
        r = sqrt(1 + b*Y) - 1;  % Short rate process
        
        % Calculate bond prices B_t = E[r_t/(e^{-rho·τ}·r_{t+τ})]
        B = zeros(1, N+1 - tau/dt);
        for k = 1:N+1 - tau/dt
            B(k) = mean(r(:,k) ./ (exp(-rho(i)*tau) * r(:,k + tau/dt)));
        end
        
        % Calculate consumption paths
        c1 = (2/b) * (sqrt(1 + b*Y) - 1);
        c2 = Y - c1;
        
        % Store results
        results{i}.Y = Y;
        results{i}.B = B;
        results{i}.c1 = c1;
        results{i}.c2 = c2;
        results{i}.t = t;
        results{i}.tau_idx = tau/dt;
    end

    %% Visualization (Adapted Style)
    tinit = 2;  % Skip initial points
    tend = min(500, N+1);  % Don't exceed array bounds
    
    % Create figure with 2x2 subplot layout
    figure('Name', 'Bond Pricing Results', 'Position', [100, 100, 1200, 800]);
    
    % Subplot 1: Sample paths of Y_t (first calibration)
    subplot(2,2,1);
    plot(results{1}.t(tinit:tend), results{1}.Y(1:5,tinit:tend));
    title('Sample Paths of Endowment (Y_t)');
    xlabel('Time');
    ylabel('Y_t');
    grid on;
    
    % Subplot 2: Bond prices comparison
    subplot(2,2,2);
    hold on;
    for i = 1:3
        plot(results{i}.t(1:end-results{i}.tau_idx), results{i}.B, 'LineWidth', 1.5);
    end
    title('Bond Prices (B_t) Across Calibrations');
    xlabel('Time');
    ylabel('B_t');
    legend(names);
    grid on;
    
    % Subplot 3: Consumption paths (mean across simulations)
    subplot(2,2,3);
    hold on;
    for i = 1:3
        plot(results{i}.t(tinit:tend), mean(results{i}.c1(:,tinit:tend)), '--', 'LineWidth', 1.5);
        plot(results{i}.t(tinit:tend), mean(results{i}.c2(:,tinit:tend)), '-', 'LineWidth', 1.5);
    end
    title('Average Consumption Paths');
    xlabel('Time');
    ylabel('Consumption');
    legend('Case 1 c1', 'Case 1 c2', 'Case 2 c1', 'Case 2 c2', 'Case 3 c1', 'Case 3 c2');
    grid on;
    
    % Subplot 4: Final distribution of Y_T
    subplot(2,2,4);
    hold on;
    for i = 1:3
        histogram(results{i}.Y(:,end), 50, 'Normalization', 'pdf');
    end
    title('Final Endowment Distribution (Y_T)');
    xlabel('Y_T');
    ylabel('Density');
    legend(names);
    grid on;
    
    %% Save Figure (Landscape PDF)
    h = gcf;
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'normalized');
    set(gcf, 'PaperPosition', [0 0 1 1]);
    print(h, '-dpdf', 'Adapted_Bond_Pricing_Results.pdf');
    
    %% Save Workspace Data
    save('Adapted_Bond_Pricing_Data.mat', 'results', 'names', 'lambda', 'rho', 'mu', 'sigma', 'gamma1');
end