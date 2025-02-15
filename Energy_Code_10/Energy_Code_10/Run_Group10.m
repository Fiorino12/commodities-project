%%                       ENERGY PROJECT - GROUP 10                       %%
%                               HJM NIG model                             %
                           
                            % Edoardo Del Bianco
                            % Mattia Fioravanti
                            % Emanuele Frigerio 
                            % Giovanni Frontali
                            % Stefano Garagiola

clear all
close all
clc
%%
addpath("Utilities\")
addpath("Calibration\")
addpath("Lookback\")
addpath("BasketOption\")
rng(10)

%% %%%%%%%%%%%%%%%%   0) PRELIMINARY: IMPORT DATA  %%%%%%%%%%%%%%%%%%%%% %%   
% Import futures prices, implied volatility surfaces, and 
% discount factors from the provided Excel file and computes interest rates and 
% option prices. The function 'ImportData' also visualizes
% the implied volatility and option price surfaces.
filename = 'DATA_FREEX.xlsx';
[strike_2026,tenor_2026,rates_2026,surf_prices2026,F0_2026,...
 strike_2028,tenor_2028,rates_2028,surf_prices2028,F0_2028,...
 rate]=ImportData(filename);

%% %%%%%%%%%%%%%%%%%%%%%%   3) CALIBRATION   %%%%%%%%%%%%%%%%%%%%%%%%%%% %%
disp('Performing calibration on 2026 surface')
% Calibrate the model on the 2026 French option prices 
% by minimizing the distance between model and market prices

% Building the system for positivity constraints on the parameters k_IG and sigma
A=[-1,0,0,0;0,-1,0,0];
b=[0;0];

% Selecting the meaningful elements of the surface
strike_2026_new=strike_2026(4:14);
surf_prices2026_new=surf_prices2026(:,4:14);

%% Minimizing the distance between model and market prices via fmincon
%!!!  Unfortunately, fmincon results in an unstable and unacceptable outcome !!!

% tic
% param_2026 = fmincon(@(param) loss_function(strike_2026_new,tenor_2026,rates_2026,surf_prices2026_new,F0_2026,param(1),...
%     param(2),param(3),param(4)),[0.05,0.9,-1.3,1],A,b,[],[],[],[],@(param) nonlincon(param));
% toc

%% Minimizing the distance between model and market prices via lsqnonlin
% ** Impose the positivity constraints on the parameters k_IG and sigma
% ** Impose the nonlinear constraint related to the drift condition Y^2sigma^2k+2thetaYk≤1
tic
param_2026 = lsqnonlin(@(param) loss_function(strike_2026_new,tenor_2026,rates_2026,surf_prices2026_new,F0_2026,param(1),...
    param(2),param(3),param(4)),[0.05,0.9,-1.3,1],[],[],A,b,[],[],@(param) nonlincon(param));
toc
%% Check Surface
% Comparison between model prices and market prices
[k_IG_2026, sigma_2026, theta_2026, Y_2026, MSE_2026] = CheckSurface(param_2026, strike_2026_new, tenor_2026, rates_2026, surf_prices2026_new, F0_2026);

%% %%%%%%%%%%%%%%%%%%%%%%%   4) SIMULATION   %%%%%%%%%%%%%%%%%%%%%%%%%%% %%
disp('Performing Simulation and pricing for Lookback Option')
% Market/Contract Input
K=300; T=1;
r=rate(15); % select the rate corresponding to the exercise date

% Montecarlo Input
M=round(252*T);             % Monitoring -> daily
Nsim=1e6;                   % Number of Simulations
dt=T/M;                     % Time step

% drift computation
phi_NIG=@(v) (1/k_IG_2026)*(1-sqrt(1+v.^2*sigma_2026^2*k_IG_2026-2*1i*theta_2026*v*k_IG_2026));
char_exp =@(v) (phi_NIG(v*Y_2026) - 1i*v*phi_NIG(-1i*Y_2026));
drift =-phi_NIG(-1i*Y_2026);

%% Simulation and pricing of the Lookback Price 
% Traditional Monte Carlo simulation without variance reduction techniques
%[price,CI]=LookbackPrice(r,F0_2026,K,dt,T,drift,k_IG_2026, sigma_2026, theta_2026, Y_2026,Nsim,M);

%% Simulation (using Antithetic Variable_S) 
% Monte Carlo simulation using antithetic variables on S
%[price_AVS,CI_AVS]=LookbackPrice_AVS(r,F0_2026,K,dt,T,drift,k_IG_2026, sigma_2026, theta_2026, Y_2026,Nsim,M);

% !!! In terms of CI_length, this technique performs slightly better than the 
% standard case but not significantly !!!
%% Simulation (using Antithetic Variable_Z)
% Monte Carlo simulation using antithetic variables on Z
[price_AVZ,CI_AVZ]=LookbackPrice_AVZ(r,F0_2026,K,dt,T,drift,k_IG_2026, sigma_2026, theta_2026, Y_2026,Nsim,M);

% !!! This is the best-performing technique  !!!

%% %%%%%%%%%%%%%%%%%%   5) CALIBRATION 2026/2028   %%%%%%%%%%%%%%%%%%%%% %%
disp('Performing joint calibration on 2026 and 2028')
% Calibrate the model on the 2026_2028 French option prices 
% by minimizing the distance between model and market prices

% Building the system for positivity constraints on the parameters k_IG and sigma
A=[-1,0,0,0;0,-1,0,0];
b=[0;0];
LB = [0, 0, -10, -10];
UB = [10, 10, 10, 10];
X0 = [0.000033100732883   4.352359857256033  -2.833097661035060   0.153132945217703]'; 

% Selecting the meaningful elements of the surface
strike_2028_new=strike_2028(1:11);
tenor_2028_new=tenor_2028(1:end-3);
surf_prices2028_new=surf_prices2028(1:end-3,1:11);

%% Minimizing the distance between model and market prices via fmincon
% ** Impose the positivity constraints on the parameters k_IG and sigma
% ** Impose the nonlinear constraint related to the drift condition Y^2sigma^2k+2thetaYk≤1
% tic
% param_26_28 = fmincon(@(param) loss_function_joint(strike_2026_new,strike_2028_new,tenor_2026,tenor_2028_new,rates_2026,rates_2028,...
%     surf_prices2026_new,surf_prices2028_new,F0_2026,F0_2028,param(1),param(2),param(3),param(4)),X0 ,A,b,[],[],LB,UB,@(param) nonlincon(param));
% toc

%% Minimizing the distance between model and market prices via sqnonlin
% ** Impose the positivity constraints on the parameters k_IG and sigma
% ** Impose the nonlinear constraint related to the drift condition Y^2sigma^2k+2thetaYk≤1
tic
param_26_28 = lsqnonlin(@(param) loss_function_joint(strike_2026_new,strike_2028_new,tenor_2026,tenor_2028_new,rates_2026,rates_2028,...
    surf_prices2026_new,surf_prices2028_new,F0_2026,F0_2028,param(1),param(2),param(3),param(4)),X0,LB,UB,A,b,[],[],@(param) nonlincon(param));
toc

%% Check Surface
% Comparison between model prices and market prices 2026
[k_IG_26_28, sigma_26_28, theta_26_28, Y_26_28, MSE_2026_conj] = CheckSurface(param_26_28, strike_2026_new, tenor_2026, rates_2026, surf_prices2026_new, F0_2026);

% Comparison between model prices and market prices 2026
[k_IG_26_28, sigma_26_28, theta_26_28, Y_26_28, MSE_2026_conj] = CheckSurface(param_26_28, strike_2028_new, tenor_2028_new, rates_2028, surf_prices2028_new, F0_2028);

%% %%%%%%%%%%%%%%%%%%%%   6) FACULTATIVE POINT   %%%%%%%%%%%%%%%%%%%%%%% %% 
disp('Performing calibration on 2028')
%% Calibration F2 (2028)
% Calibrate the model on the 2028 French option prices 
% by minimizing the distance between model and market prices

% Building the system for positivity constraints on the parameters k_IG and sigma
A=[-1,0,0,0;0,-1,0,0];
b=[0;0];

%% Minimizing the distance between model and market prices via fmincon
%!!!  Unfortunately, fmincon results in an unstable and unacceptable outcome !!!

%tic
%param_2028 = fmincon(@(param) loss_function(strike_2028_new,tenor_2028_new,rates_2028,surf_prices2028_new,F0_2028,param(1),...
%    param(2),param(3),param(4)),[0.05,0.9,-1.3,1],A,b,[],[],[],[],@(param) nonlincon(param));
%toc

%% Minimizing the distance between model and market prices via lsqnonlin
% ** Impose the positivity constraints on the parameters k_IG and sigma
% ** Impose the nonlinear constraint related to the drift condition Y^2sigma^2k+2thetaYk≤1

tic
param_2028 = lsqnonlin(@(param) loss_function(strike_2028_new,tenor_2028_new,rates_2028,surf_prices2028_new,F0_2028,param(1),...
    param(2),param(3),param(4)),[0.05,0.9,-1.3,1],[],[],A,b,[],[],@(param) nonlincon(param));
toc

%% Check Surface 2028
% Comparison between model prices and market prices
[k_IG_2028, sigma_2028, theta_2028, Y_2028, MSE_2028] = CheckSurface(param_2028, strike_2028_new, tenor_2028_new, rates_2028, surf_prices2028_new, F0_2028);

%% Option Pricing
% disp('Simulation and Pricing Best of Option')
% drift computation
drift_2026=drift;

phi_NIG=@(v) (1/k_IG_2028)*(1-sqrt(1+v.^2*sigma_2028^2*k_IG_2028-2*1i*theta_2028*v*k_IG_2028));
char_exp=@(v) (phi_NIG(v*Y_2028) - 1i*v*phi_NIG(-1i*Y_2028));
drift_2028=-phi_NIG(-1i*Y_2028);

%% Basket Option Pricing via MC simulation 
%[price_BO,CI_BO]=Pricing_BasketOption(r,F0_2026,K,T,drift_2026,k_IG_2026, sigma_2026, theta_2026, Y_2026,...
%    drift_2028,F0_2028,k_IG_2028, sigma_2028, theta_2028, Y_2028,Nsim);

%% Basket Option Pricing via MC simulation (using Antithetic Variable_S)
%[price_BO_AVS,CI_BO_AVS]=Pricing_BasketOption_AVS(r,F0_2026,K,T,drift_2026,k_IG_2026, sigma_2026, theta_2026, Y_2026,...
%    drift_2028,F0_2028,k_IG_2028, sigma_2028, theta_2028, Y_2028,Nsim);

%% Basket Option Pricing via MC simulation (using Antithetic Variable_Z)
[price_BO_AVZ,CI_BO_AVZ]=Pricing_BasketOption_AVZ(r,F0_2026,K,T,drift_2026,k_IG_2026, sigma_2026, theta_2026, Y_2026,...
    drift_2028,F0_2028,k_IG_2028, sigma_2028, theta_2028, Y_2028,Nsim);

%% Basket under the revised framework, point 5
ClosedPrice = Pricing_BasketOption_ClosedFormula(r, K, F0_2026, F0_2028 , T, ...
    param_26_28(1), param_26_28(2), param_26_28(3), param_26_28(4)); 
