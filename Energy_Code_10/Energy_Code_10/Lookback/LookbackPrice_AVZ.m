function [price_AVZ,CI_AVZ]=LookbackPrice_AVZ(r,F0,K,dt,T,drift,k_IG, sigma, theta, Y,Nsim,M)

% LookbackPrice_AVS computes the price of a lookback option using AntitheticVariable.
% The pricing is done via Monte Carlo simulations.

% Input:
%   r        : Risk-free interest rate (exercise date)
%   F0       : Initial price of the underlying asset
%   K        : Strike price of the option
%   dt       : Time step for the simulation
%   T        : Maturity of the option (time to expiry)
%   drift    : Drift component associated to the integral on A
%   k_IG,sigma, theta,Y : NIG params
%   Nsim     : Number of Monte Carlo simulations
%   M        : Monitoring

% Output:
%   price_AVS : Estimated price of the lookback option
%   CI_AVS    : Confidence interval for the option price

    %% Simulation 
    % Vector Inizialization
    % Variable1
    X_old_1 = zeros(Nsim, 1);
    X_new_1 = zeros(Nsim, 1);
    X_max_1 = zeros(Nsim, 1);

    % Variable2
    X_old_2 = zeros(Nsim, 1);
    X_new_2 = zeros(Nsim, 1);
    X_max_2 = zeros(Nsim, 1);

    Z = randn(Nsim, M);  % Simulation Standard Gaussian

    % Monte Carlo Simulation Loop
    for j = 1:M
        [dS, ~] = fastInverseGaussian(dt, dt^2 / k_IG, Nsim);
    
        % Antithetic variable Z standard Gaussian:
        X_new_1 = X_old_1 + drift * dt + ...
            Y * (theta * dS + sigma * sqrt(dS) .* Z(:, j));
        X_max_1 = max(X_new_1, X_max_1);
        X_old_1 = X_new_1;

        X_new_2 = X_old_2 + drift * dt + ...
            Y * (theta * dS - sigma * sqrt(dS) .* Z(:, j));
        X_max_2 = max(X_new_2, X_max_2);
        X_old_2 = X_new_2;
    end

    % Check Variable1
    disp('Martingality check first variable: must be one')
    fprintf('E[exp(X1t)]: %d\n', mean(exp(X_old_1)))

    % Check Variable2
    disp('Martingality check second variable: must be one')
    fprintf('E[exp(X2t)]: %d\n', mean(exp(X_old_2)))

    % Pricing   
    Payoff_disc_AVZ_1 = max(F0 * exp(X_max_1) - K, 0) * exp(-r * T);
    Payoff_disc_AVZ_2 = max(F0 * exp(X_max_2) - K, 0) * exp(-r * T);
    [price_AVZ, ~, CI_AVZ] = normfit((Payoff_disc_AVZ_1 + Payoff_disc_AVZ_2) / 2);

    % Create table for Price and Confidence Interval (CI)
    results_table = table(price_AVZ, CI_AVZ(1), CI_AVZ(2), CI_AVZ(2)-CI_AVZ(1) , 'VariableNames', {'Price_AVZ', 'CI_Lower', 'CI_Upper', 'CI_length'});

    % Display the results table
    disp('Lookback Option Pricing Results:')
    disp(results_table)

    %% check call option
    % 2026
    payoff_call_2026 = (max(F0 * exp(X_old_1)-K,0)*exp(-r*T)+max(F0 * exp(X_old_2)-K,0)*exp(-r*T))/2;
    [price_call_2026_MC,~,CI_call_2026]=normfit(payoff_call_2026);

    phi_NIG=@(v) T*(1/k_IG)*(1-sqrt(1+v.^2*sigma^2*k_IG-2*1i*theta*v*k_IG));
    CharFunc=@(v) exp(phi_NIG(v*Y) - 1i*v*phi_NIG(-1i*Y));
    price_call_2026_FFT = FFT_CM_Call_NIG(K,T,r,F0,CharFunc);

    call_table_2026 = table(price_call_2026_MC,price_call_2026_FFT,'VariableNames', {'MC 2026', 'Carr-Madan 2026'});
    disp('Call Pricing Results:')
    disp(call_table_2026)


end