function [price_6,CI_6]=Pricing_BasketOption_AVZ(r,F0_2026,K,T,drift_2026,k_IG_2026, sigma_2026, theta_2026, Y_2026,...
    drift_2028,F0_2028,k_IG_2028, sigma_2028, theta_2028, Y_2028,Nsim)

% Pricing_ex6 Computes the price of an option with payoff
% max(max(F2026(1y),F2028(1y))-K,0)
%
% INPUTS:
%   r          : Risk-free interest rate (exercise date)
%   F0_2026    : Price at time 0 for the 2026 swap.
%   K          : Strike price of the option.
%   T          : Maturity time of the option.
%   drift_2026, k_IG_2026, sigma_2026, theta_2026, Y_2026 : process parameters (2026 swap)
%   drift_2028, k_IG_2028, sigma_2028, theta_2028, Y_2028 : process parameters (2028 swap)
%   Nsim       - Number of Monte Carlo simulations.
%
% OUTPUTS:
%   price_6    - Estimated option price.
%   CI_6       - Confidence interval for the estimated price.

    %% Simulation swap2026
    Z_2026=randn(Nsim,1);

    [dS_2026,~] = fastInverseGaussian(T, T^2/k_IG_2026, Nsim);
    XT_2026_1=drift_2026*T+Y_2026*(theta_2026*dS_2026+sigma_2026*sqrt(dS_2026).*Z_2026);
    XT_2026_2=drift_2026*T+Y_2026*(theta_2026*dS_2026-sigma_2026*sqrt(dS_2026).*Z_2026);
    FT_2026_1=F0_2026*exp(XT_2026_1);
    FT_2026_2=F0_2026*exp(XT_2026_2);
    % Martingality check swap 2026
    disp('Martingality check swap 2026: must be one')
    fprintf('E[exp(Xt_1)]: %d\n', mean(exp(XT_2026_1)))
    fprintf('E[exp(Xt_2)]: %d\n', mean(exp(XT_2026_2)))
    
    %% Simulation swap2028
    Z_2028=randn(Nsim,1);

    [dS_2028,~] = fastInverseGaussian(T, T^2/k_IG_2028, Nsim);
    XT_2028_1=drift_2028*T+Y_2028*(theta_2028*dS_2028+sigma_2028*sqrt(dS_2028).*Z_2028);
    XT_2028_2=drift_2028*T+Y_2028*(theta_2028*dS_2028-sigma_2028*sqrt(dS_2028).*Z_2028);
    FT_2028_1=F0_2028*exp(XT_2028_1);
    FT_2028_2=F0_2028*exp(XT_2028_2);
    % Martingality check swap 2028
    disp('Martingality check swap 2026: must be one')
    fprintf('E[exp(Xt_1)]: %d\n', mean(exp(XT_2028_1)))
    fprintf('E[exp(Xt_2)]: %d\n', mean(exp(XT_2028_2)))
    
    %% Option Pricing
    discpayoff_1=max(max(FT_2026_1,FT_2028_1)-K,0)*exp(-r*T);
    discpayoff_2=max(max(FT_2026_2,FT_2028_2)-K,0)*exp(-r*T);
    [price_6,~,CI_6]=normfit((discpayoff_1+discpayoff_2)/2);

    % Create table for Price and Confidence Interval (CI)
    results_table = table(price_6, CI_6(1), CI_6(2),CI_6(2)-CI_6(1), 'VariableNames', {'Price', 'CI_Lower', 'CI_Upper', 'CI_length'});

    % Display the results table
    disp('Option Pricing Results:')
    disp(results_table)

    %% check call option
    % 2026
    payoff_call_2026 = (max(FT_2026_1-K,0)*exp(-r*T)+max(FT_2026_2-K,0)*exp(-r*T))/2;
    [price_call_2026_MC,~,CI_call_2026]=normfit(payoff_call_2026);

    phi_NIG=@(v) T*(1/k_IG_2026)*(1-sqrt(1+v.^2*sigma_2026^2*k_IG_2026-2*1i*theta_2026*v*k_IG_2026));
    CharFunc=@(v) exp(phi_NIG(v*Y_2026) - 1i*v*phi_NIG(-1i*Y_2026));
    price_call_2026_FFT = FFT_CM_Call_NIG(K,T,r,F0_2026,CharFunc);

    call_table_2026 = table(price_call_2026_MC,price_call_2026_FFT,'VariableNames', {'MC 2026', 'Carr-Madan 2026'});

    % 2028
    payoff_call_2028 = (max(FT_2028_1-K,0)*exp(-r*T)+max(FT_2028_2-K,0)*exp(-r*T))/2;
    [price_call_2028_MC,~,CI_call_2028]=normfit(payoff_call_2028);

    phi_NIG=@(v) T*(1/k_IG_2028)*(1-sqrt(1+v.^2*sigma_2028^2*k_IG_2028-2*1i*theta_2028*v*k_IG_2028));
    CharFunc=@(v) exp(phi_NIG(v*Y_2028) - 1i*v*phi_NIG(-1i*Y_2028));
    price_call_2028_FFT = FFT_CM_Call_NIG(K,T,r,F0_2028,CharFunc);

    call_table_2028 = table(price_call_2028_MC,price_call_2028_FFT,'VariableNames', {'MC 2028', 'Carr-Madan 2028'});

    disp('Call Pricing Results:')
    disp(call_table_2026)
    disp(call_table_2028)
  
end