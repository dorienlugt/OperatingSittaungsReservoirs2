%%  Reservoir details

function [ResDetails] = reservoirDetails(optimization)

    % Number of the reservoir(s) included
    Index = optimization.Index;
    % Length of one timestep (s)
    dt = optimization.dt;       
    
%%  Reservoir characteristics
    Reservoir_names = string({'Yenwe' 'Phyu' 'Paunglaung' 'Kabaung' ...
        'Baingda' 'Sinthe' 'Kawliya' 'Bawni' 'Swa' 'Yanaungmyin' ...
        'Chaungmange' 'Ngalaik' 'Madan' 'Yezin' 'Pathi' 'Yathoe' ...
        'Chaungmagyi_pyi' 'Swegyin' 'Thaukyegat' 'Minye' 'Kunchaung'}); 
    
    % Reservoir catchment areas (m2) 
    Res_Catch_area_All = 10^6*[804 1051 4791 1178 246 809 116 45 1077 ...
        28 265 328 100 33 66 33 119 875 2161 23 875];     
    Res_Catch_area = Res_Catch_area_All;
    % Conduit capacity of the reservoirs (m3 /s)
    Reservoir_Conduit = dt*[50 50 306 39 33.98 17 25 4.25 21.24 1.81 ...
        7.08 8.67 7.93 8.39 3.54 2.83 1.12 213 210 0.19 66];      
    % Storage capacity of the reservoirs (m3)
    Reservoir_Cap = 10^6*[988 726 360 890 461 176 174 38 267 19 113 ...
        93 45 90 38 19 12 1450 295 2 843]; 
    % Irrigable area of the reservoirs (m2)
    Res_Irr_Areas = 10^6*[480 405 142 216 189 131 34 13 142 8 ...
        32 85 32 65 12 8 12 254 627 4 254];                               
    % Waterspread area of the reservoirs (m2)
    Res_Waterspread = 10^4*[7433 5896 1660 5931 623 2700 1621 782 ...
        2792 665 1274 1558 450 1125 376 567 1673 9503 2367 79 5753];  
    % Number of subcatchment the reservoirs is located in
    Reservoir_Catchm = [16 9 1 9 17 1 17 17 3 2 2 2 2 1 6 9 2 15 7 9 12];
    % Area of the subcatchments (m2)
    Catchm_areas = 10^6*[3014+4920 248+3261 1517 240 295 344 2473 326 ...
        3421 962 491 1903 399 642 1732 1433 2450];               

    %   Select details of included reservoirs 
    if (length(Index)<length(Reservoir_Catchm))
        Res_Catch_area = Res_Catch_area_All(Index);
        Reservoir_Catchm = Reservoir_Catchm(Index);
        Reservoir_Conduit = Reservoir_Conduit(Index);
        Reservoir_Cap = Reservoir_Cap(Index);
        Res_Irr_Areas = Res_Irr_Areas(Index);
        Res_Waterspread = Res_Waterspread(Index);       
        Catchm_areas = Catchm_areas(Reservoir_Catchm);
    end  
    
    ResDetails.Names = Reservoir_names;
    ResDetails.S = Reservoir_Cap';                  
    ResDetails.G = Reservoir_Conduit';                    
    ResDetails.ResA = Res_Waterspread;
    ResDetails.CatchmA = Catchm_areas;
    ResDetails.ResCatchmA = Res_Catch_area;
    ResDetails.ResCatchmA_All = Res_Catch_area_All;
    ResDetails.ResCatchm = Reservoir_Catchm;
  
%%  Irrigation demands
    % Water demands in (m/day) for Rice in Myanmar for Jan-Dec
    ETcrop = 10^(-3)*[3.1 3.1 2.4 2.1 0 3.2 3.5 1.3 1.6 1.6 1.8 4.9];
    % Irrigation demands (m3/day) per reservoir (columns) per month (rows) 
    ResDetails.Irr = ETcrop'*Res_Irr_Areas;                            
    
%%  Hydropower parameters
    % Hydropower constant factor to calculate
    % Energy (MWh) = HPc * Conduit outflow * Water level reservoir 
    HPc = 0.4*9.81*10^(-3)*(1/3600);   
    ResDetails.HPc = HPc;
   
end