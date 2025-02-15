function [k_IG, sigma, theta, Y, MSE] = CheckSurface(param, strike, tenor, rates, surf_prices, F0)

% This function computes the model price surface for options based on the
% Normal Inverse Gaussian (NIG) process using parameters provided in the 
% input. It calculates the Mean Squared Error (MSE) between the model prices 
% and the market prices and plots the two surfaces for comparison.
%
% Input:
%   param       : NIG model parameters:
%                   param(1) = k_IG
%                   param(2) = sigma
%                   param(3) = theta
%                   param(4) = Y
%   strike      : Strike prices for the options.
%   tenor       : Option tenors.
%   rates       : Interest rates corresponding to the tenors.
%   surf_prices : A matrix of market option prices.
%   F0          : The current forward price of the underlying asset.
%
% Output:
%   k_IG,sigma,theta,Y : NIG parameters.    
%   MSE                : The Mean Squared Error between model prices and market prices.

    %% Compute MSE and model prices
    % Extract parameters from the input vector
    k_IG = param(1);
    sigma = param(2);
    theta = param(3);
    Y = param(4);
    
    % Create a table to display the parameters
    param_table = table(k_IG, sigma, theta, Y, ...
        'VariableNames', {'k_IG', 'sigma', 'theta', 'Y'});
    disp('NIG parameters');
    disp(param_table);


    % Compute the Mean Squared Error (MSE) using the loss function
    MSE = loss_function(strike, tenor, rates, surf_prices, F0, k_IG, sigma, theta, Y);
    % Print MSE
    fprintf('Mean Squared Error %d \n',MSE);

    % Initialize the model price surface
    model_price = zeros(length(tenor), length(strike));

    % Loop through each tenor to calculate model prices
    for jj = 1:length(tenor)
        % Define the characteristic function components
        phi_NIG = @(v) tenor(jj) * (1 / k_IG) * (1 - sqrt(1 + v.^2 * sigma^2 * k_IG - 2 * 1i * theta * v * k_IG));
        CharFunc = @(v) exp(phi_NIG(v * Y) - 1i * v * phi_NIG(-1i * Y));

        % Compute model prices using FFT-based method
        model_price(jj, :) = FFT_CM_Call_NIG(strike, tenor(jj), rates(jj), F0, @(v) CharFunc(v));
    end

    %% Plot the model and market price surfaces
    figure()
    surf(strike, tenor, model_price, 'FaceColor', 'r')
    hold on
    surf(strike, tenor, surf_prices, 'FaceColor', 'b')
    title('Comparison of Model and Market Prices', 'FontSize', 14, 'FontWeight', 'bold')
    xlabel('Strike Price', 'FontSize', 12)
    ylabel('Tenor (Years)', 'FontSize', 12)
    zlabel('Option Price', 'FontSize', 12)
    legend("Model", "Market")
    grid on
end