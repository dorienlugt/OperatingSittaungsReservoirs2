%%  Load data 
%   External disturbances are outflows from subcatchments and 
%   reservoir catchments
%   Based on inflow dataseries Paunglaung 

function [Qm,Qseas,Day,Month,Year,Area] = ...
    forcingdata_Seasonality(optimization,ResDetails)

    T_start = optimization.T_start;
    T_opt = optimization.T_opt;
    T_shift = optimization.T_shift; 
    T_mov = optimization.T_mov;

    Reservoir_names = ResDetails.Names;      
    ResCatchmA_All = ResDetails.ResCatchmA_All;

    Area = ResCatchmA_All(Reservoir_names == 'Paunglaung');
    
    Qseas = csvread('Input\Paunglaung3_Seasonality.csv',0,4);    
    Qm = csvread('Input\Paunglaung3.csv',0,5);
    Qm(isnan(Qm)) = 0;
    Month = csvread('Input\Paunglaung3.csv',0,1);
    Day = csvread('Input\Paunglaung3.csv',0,0);
    Year = csvread('Input\Paunglaung3.csv',0,2);
 
    Qm = Qm((T_start-1+(1:(T_opt+(T_mov+1)*T_shift))),1);
    Qseas = Qseas((T_start-1+(1:(T_opt+(T_mov+1)*T_shift))),1);
    Month = Month((T_start-1+(1:(T_opt+(T_mov+1)*T_shift))),1);
    Year = Year((T_start-1+(1:(T_opt+(T_mov+1)*T_shift))),1);
    Day = Day((T_start-1+(1:(T_opt+(T_mov+1)*T_shift))),1);
    
end