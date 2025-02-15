function [price,CI]=LookbackPrice(r,F0,K,dt,T,drift,k_IG, sigma, theta, Y,Nsim,M)

% LookbackPrice computes the price of a lookback option.
% The pricing is done via Monte Carlo simulations.
% At the end we compute also the price for a Call option using MC simulation and Carr-Madan Algorithm to check our
% simulation 

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
%   price     : Estimated price of the lookback option
%   CI        : Confidence interval for the option price

    %% Simulation 
    % Vector Inizialization
    % Variable1
    X_old = zeros(Nsim, 1);
    X_new = zeros(Nsim, 1);
    X_max = zeros(Nsim, 1);

    Z = randn(Nsim, M);  % Simulation Standard Gaussian

    % Monte Carlo Simulation Loop
    for j = 1:M
        % IG Simulation 
        [dS, ~]= fastInverseGaussian(dt, dt^2 / k_IG, Nsim);
        X_new = X_old + drift * dt + ...
            Y * (theta * dS + sigma * sqrt(dS) .* Z(:, j));
        X_max = max(X_new, X_max);
        X_old = X_new;

    end

    % Martingality check
    disp('Martingality check first variable: must be one')
    fprintf('E[exp(Xt)]: %d\n', mean(exp(X_old)))

    % Pricing   
    Payoff_disc = max(F0 * exp(X_max) - K, 0) * exp(-r * T);
    [price, ~, CI] = normfit(Payoff_disc);

    % Create table for Price and Confidence Interval (CI)
    results_table = table(price, CI(1), CI(2),CI(2)-CI(1), 'VariableNames', {'Price', 'CI_Lower', 'CI_Upper', 'CI_length'});

    % Display the results table
    disp('Lookback Option Pricing Results:')
    disp(results_table)

    %% check call option
    % call price using MC simulation
    discpayoff_call=exp(-r*T) * max(F0 * exp(X_old) - K, 0);
    [price_call_MC,~,CI_call_MC]=normfit(discpayoff_call);
    
    % call price using Carr-Madan algorithm
    phi_NIG=@(v) T*(1/k_IG)*(1-sqrt(1+v.^2*sigma^2*k_IG-2*1i*theta*v*k_IG));
    CharFunc=@(v) exp(phi_NIG(v*Y) - 1i*v*phi_NIG(-1i*Y));
    
    price_call_CM=FFT_CM_Call_NIG(K,T,r,F0,@(v) CharFunc(v));

    % Create table for Price and Confidence Interval (CI)
    results_table_call = table(price_call_MC, price_call_CM, 'VariableNames', {'Call price with Monte Carlo', 'Call price with Carr-Madan'});

    % Display the results table
    disp('Call Option Pricing Results:')
    disp(results_table_call)

end