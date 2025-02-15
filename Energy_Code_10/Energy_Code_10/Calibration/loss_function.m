function mse = loss_function(strike,tenor,rates,market_price,F0,k_IG,sigma,theta,Y)

% This function calculates the mean squared error (MSE) between the model-generated option prices
% and the observed market prices. It is used to calibrate the model parameters by minimizing the loss.
%
% Inputs:
%   - strike        : A vector of strike prices.
%   - tenor         : A vector of time to maturities (tenors).
%   - rates         : A vector of interest rates corresponding to the tenors.
%   - market_price  : A matrix of observed market prices (rows correspond to tenors, columns to strikes).
%   - F0            : The initial swap price of the underlying asset.
%   - k_IG, sigma, theta, Y : NIG parameters 
%
% Output:
%   - mse           : The mean squared error (MSE) between the model prices and market prices.

mse = 0;
model_price = zeros(length(tenor),length(strike));
w = 1./(abs(strike-F0)).^2;
for jj=1:length(tenor)
    
    phi_NIG=@(v,k_IG,sigma,theta) tenor(jj)*(1/k_IG)*(1-sqrt(1+v.^2*sigma^2*k_IG-2*1i*theta*v*k_IG));
    CharFunc=@(v,k_IG,sigma,theta,Y) exp(phi_NIG(v*Y,k_IG,sigma,theta) - 1i*v*phi_NIG(-1i*Y,k_IG,sigma,theta));
    model_price(jj,:) = FFT_CM_Call_NIG(strike,tenor(jj),rates(jj),F0,@(v) CharFunc(v,k_IG,sigma,theta,Y));

    % MSE computation
    mse = mse + sum(w.*((model_price(jj,:) - market_price(jj,:))).^2);
end
mse = mse/(length(tenor)*length(strike));
end