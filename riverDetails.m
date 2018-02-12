%%  River details

function [RivDetails] = riverDetails(optimization)

    % Length of the river segments (m)
    Stream_lengths = [28273 1817+39182 8491 43068 10599 41983 ...
        8631 11286 33336 9695 23533 26010 28673 27176 27295 34484 9999];     
    % Conversion from storage to depth river segment (width = 200 m)
    Conv = 1./(Stream_lengths'.*200); 

    RivDetails.Conv = Conv;
    dt = optimization.dt;

    % Travel speed wave (m/s)
    T = 1.25;    
    % Storage-time coefficient K(s) (Muskingum routing parameter)
    K = Stream_lengths./T;         
    % Weighting factor inflow versus outflow (Muskingum routing parameter)
    X = 0.3;
    RivDetails.c1 = (dt+2*K*X)./(dt+2*K-2*K*X);
    RivDetails.c2 = (dt-2*K*X)./(dt+2*K-2*K*X);
    RivDetails.c3 = (-dt+2*K-2*K*X)./(dt+2*K-2*K*X);

end

