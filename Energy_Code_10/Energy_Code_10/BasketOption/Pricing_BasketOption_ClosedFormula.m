function ClosedPrice = Pricing_BasketOption_ClosedFormula(r, Strike, F0_2026, F0_2028 ,T, k_IG, sigma, theta, Y)
% Pricing_ex6 Computes the price of an option with payoff
% max(max(F2026(1y),F2028(1y))-K,0)
% under the point 5 framework, i.e. both swaps drove by the same lévy process 

% INPUTS:
%   r          : Risk-free interest rate (exercise date)
%   F0_2026    : Price at time 0 for the 2026 swap.
%   F0_2028    : Price at time 0 for the 2028 swap.
%   K          : Strike price of the option.
%   T          : Maturity time of the option.
%   drift, k_IG, sigma, theta, Y: process parameters under the new framwork

%
% OUTPUTS:
%   price    - Price obtained via the Carr-Madan formula

% Characteristic exponent initialization 
phi_NIG=@(v,k_IG,sigma,theta) T*(1/k_IG)*(1-sqrt(1+v.^2*sigma^2*k_IG-2*1i*theta*v*k_IG));
CharFunc=@(v,k_IG,sigma,theta,Y) exp(phi_NIG(v*Y,k_IG,sigma,theta) - 1i*v*phi_NIG(-1i*Y,k_IG,sigma,theta));

% Closed Formula computation 
ClosedPrice = FFT_CM_Call_NIG(Strike, T, r, max(F0_2026, F0_2028), @(v)CharFunc(v,k_IG,sigma,theta,Y)); 
call_table = table(ClosedPrice, 'VariableNames', { 'Carr-Madan'});

% Display the results table
disp('Option Pricing Results:')
disp(call_table)
end