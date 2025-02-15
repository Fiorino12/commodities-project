function [strike_2026,tenor_2026,rates_2026,surf_prices2026,F0_2026,...
          strike_2028,tenor_2028,rates_2028,surf_prices2028,F0_2028, rate]=ImportData(filename)

% This function imports futures prices, implied volatility surfaces, and 
% discount factors from the provided Excel file and computes interest rates and 
% option prices based. The function also visualizes the implied volatility 
% and option price surfaces.

% Input:
%   filename : the name of the Excel file containing the data. 
%              Sheet 1: Futures prices with columns 'Expiry' and 'Price'.
%              Sheet 2: Implied Volatility Surface for 2026.
%              Sheet 3: Implied Volatility Surface for 2028.
%              Sheet 4: Discount factors.
%
% Output:
%   strike_2026     : Array of strike prices for the 2026 options.
%   tenor_2026      : Array of tenors for the 2026 options.
%   rates_2026      : Interest rates corresponding to the tenors for 2026.
%   surf_prices2026 : Matrix of option prices for the 2026 surface.
%   strike_2028     : Array of strike prices for the 2028 options.
%   tenor_2028      : Array of tenors for the 2028 options.
%   rates_2028      : Interest rates corresponding to the tenors for 2028.
%   surf_prices2028 : Matrix of prices for the 2028 surface.  
%   rate            : rates curve imported from the Excel File.

    %% Import Data
    % Select Excel sheets
    fprintf('Importing data from %s file \n',filename);
    sheets = sheetnames(filename);

    % Import Futures' Prices from 1st sheet
    futures = readtable(filename, 'Sheet', sheets{1});
    expiry_futures = futures.Expiry;
    price_futures = futures.Price;

    % Import ImpliedVolatility Surface 2026 from 2nd sheet
    options_2026 = readtable(filename, 'Sheet', sheets{2});
    strike_2026 = options_2026(1,2:end);
    strike_2026 = table2array(strike_2026);                   % Strike Price
    tenor_2026 = options_2026(2:end,1);
    tenor_2026 = table2array(tenor_2026);                     % Tenor
    implied_vol_2026 = options_2026(2:end,2:end);
    implied_vol_2026 = table2array(implied_vol_2026);         % Implied Volatility

    % Import ImpliedVolatility Surface 2028 from 3rd sheet
    options_2028 = readtable(filename, 'Sheet', sheets{3});
    strike_2028 = options_2028(1,2:end);
    strike_2028 = table2array(strike_2028);                   % Strike Price
    tenor_2028 = options_2028(2:end,1);
    tenor_2028 = table2array(tenor_2028);                     % Tenor
    implied_vol_2028 = options_2028(2:end,2:end);
    implied_vol_2028 = table2array(implied_vol_2028);         % Implied Volatility

    % Import Discounts from from 4th sheet
    discounts = table2array(readtable(filename, 'Sheet', sheets{4}));                 % Discounts
    % Adjusting the discrepancy in Excel-MATLAB date conversion
    adj=datetime(discounts(1,1), 'ConvertFrom', 'datenum') - datetime('14-Nov-2024'); 
    dates = datetime(discounts(1,:), 'ConvertFrom','datenum')-adj;                    % Date

    %% Compute rates associated to option tenors
    % Set the current date
    init_date = datetime('04-Nov-2024');

    % Apply LinearInterpolation for the computation of the interest rates related
    % to each tenor
    time_interval = yearfrac(init_date, dates,3);
    rate = - log(discounts(2,:))./time_interval;
    rates_2026 = interp1(time_interval, rate,tenor_2026);
    rates_2028 = interp1(time_interval, rate,tenor_2028);

    %% Convert ImpliedVolatility Surface to PriceOption Surface
    % Plot ImpliedVolatility Surface for 2026
    figure()
    surf(strike_2026, tenor_2026, implied_vol_2026)
    title('Implied Volatility Surface (2026)', 'FontSize', 14, 'FontWeight', 'bold')
    xlabel('Strike Price', 'FontSize', 12)
    ylabel('Tenor (Years)', 'FontSize', 12)
    zlabel('Implied Volatility', 'FontSize', 12)
    grid on
    
    % Plot ImpliedVolatility Surface for 2028
    figure()
    surf(strike_2028, tenor_2028, implied_vol_2028)
    title('Implied Volatility Surface (2028)', 'FontSize', 14, 'FontWeight', 'bold')
    xlabel('Strike Price', 'FontSize', 12)
    ylabel('Tenor (Years)', 'FontSize', 12)
    zlabel('Implied Volatility', 'FontSize', 12)
    grid on

    % Convertion ImpliedVolatility -> PriceOption
    surf_prices2026=zeros(size(implied_vol_2026));
    surf_prices2028=zeros(size(implied_vol_2028));
    F0_2026 = price_futures(16,1);
    for ii=1:size(implied_vol_2026,1)
        for j=1:size(implied_vol_2026,2)
            [surf_prices2026(ii,j),~] = blkprice(F0_2026,strike_2026(j),rates_2026(ii),tenor_2026(ii), implied_vol_2026(ii,j));
        end
    end
    F0_2028 = price_futures(18,1);
    for ii=1:size(implied_vol_2028,1)
        for j=1:size(implied_vol_2028,2)
            [surf_prices2028(ii,j),~] = blkprice(F0_2028,strike_2028(j),rates_2028(ii),tenor_2028(ii), implied_vol_2028(ii,j));
        end
    end
    
    %% Plot PriceOption Surfaces
    % Plot Price Surface for 2026
    figure()
    surf(strike_2026, tenor_2026, surf_prices2026)
    title('Option Price Surface (2026)', 'FontSize', 14, 'FontWeight', 'bold')
    xlabel('Strike Price', 'FontSize', 12)
    ylabel('Tenor (Years)', 'FontSize', 12)
    zlabel('Option Price', 'FontSize', 12)
    grid on

    % Plot Price Surface for 2028
    figure()
    surf(strike_2028, tenor_2028, surf_prices2028)
    title('Option Price Surface (2028)', 'FontSize', 14, 'FontWeight', 'bold')
    xlabel('Strike Price', 'FontSize', 12)
    ylabel('Tenor (Years)', 'FontSize', 12)
    zlabel('Option Price', 'FontSize', 12)
    grid on

end