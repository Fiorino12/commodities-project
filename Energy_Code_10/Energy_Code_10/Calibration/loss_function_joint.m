function mse=loss_function_joint(strike1,strike2,tenor1,tenor2,rates1,rates2,market_price1,market_price2,F0_1,F0_2,k_IG,sigma,theta,Y)

% This function calculates the mean squared error (MSE) between the model-generated option prices
% and the observed market prices (2026-2028 surfaces). It is used to calibrate the model parameters 
% by minimizing the loss.


% Inputs:
%   - strike1(2)          : A vector of strike prices.
%   - tenor1(2)           : A vector of time to maturities (tenors).
%   - rates1(2)           : A vector of interest rates corresponding to the tenors.
%   - market_price1(2)    : A matrix of observed market prices (rows correspond to tenors, columns to strikes).
%   - F0_1(_2)            : The initial swap price of the underlying asset.
%   - k_IG, sigma, theta, Y : NIG parameters 
%
% Output:
%   - mse           : The mean squared error (MSE) between the model prices and market prices.

    mse = 0;
    model_price1 = zeros(length(tenor1),length(strike1));
    w1 = 1./(abs(strike1-F0_1)).^2;
    for jj=1:length(tenor1)
        phi_NIG=@(v,k_IG,sigma,theta) tenor1(jj)*(1/k_IG)*(1-sqrt(1+v.^2*sigma^2*k_IG-2*1i*theta*v*k_IG));
        CharFunc=@(v,k_IG,sigma,theta,Y) exp(phi_NIG(v*Y,k_IG,sigma,theta) - 1i*v*phi_NIG(-1i*Y,k_IG,sigma,theta));
        model_price1(jj,:) =FFT_CM_Call_NIG(strike1,tenor1(jj),rates1(jj),F0_1,@(v) CharFunc(v,k_IG,sigma,theta,Y));
        mse = mse + sum(w1.*((model_price1(jj,:) - market_price1(jj,:))).^2);
    end
    model_price2 = zeros(length(tenor2),length(strike2));
    w2 = 1./(abs(strike2-F0_2)).^2;
    for jj=1:length(tenor2)
        phi_NIG=@(v,k_IG,sigma,theta) tenor2(jj)*(1/k_IG)*(1-sqrt(1+v.^2*sigma^2*k_IG-2*1i*theta*v*k_IG));
        CharFunc=@(v,k_IG,sigma,theta,Y) exp(phi_NIG(v*Y,k_IG,sigma,theta) - 1i*v*phi_NIG(-1i*Y,k_IG,sigma,theta));
        model_price2(jj,:) =FFT_CM_Call_NIG(strike2,tenor2(jj),rates2(jj),F0_2,@(v) CharFunc(v,k_IG,sigma,theta,Y));
        mse = mse + sum(w2.*((model_price2(jj,:) - market_price2(jj,:))).^2);
    end
    mse = mse/((length(tenor1)+length(tenor2))*(length(strike1)+length(strike2)));
end