%   Bounds on the variables and constraints

function [lb,ub,cl,cu] = bounds(auxdata)
    
    % Setting for the optimization
    T_opt = auxdata.optimization.T_opt;
    K_Flood = auxdata.optimization.flood.K_Flood;           
	lev_Flood = auxdata.optimization.flood.lev_Flood;	
	seg_Flood = auxdata.optimization.flood.seg_Flood;
    consIncluded = auxdata.optimization.consIncluded;      
    
    % Forcing data
    Qm_opt = auxdata.forcing.Qm_opt;
    Month_opt = auxdata.forcing.Month_opt;
    Area = auxdata.forcing.Area;
    
    % Reservoir details     
    S = auxdata.ResDetails.S;      
    G = auxdata.ResDetails.G;  
    Irr = auxdata.ResDetails.Irr;                    
    CatchmA = auxdata.ResDetails.CatchmA;
    
    % River details 
	Conv = auxdata.RivDetails.Conv;
    
    % System details
    J = auxdata.system.J;
    I = auxdata.system.I;  
    D = auxdata.system.D;
    K = auxdata.system.K;
    
%%  Variable bounds: lb < x < ub
    lb = [zeros(J,1); -Inf*zeros(3*I,1); zeros(3*J,1); zeros(K_Flood,1)];          
    ub = [S; Inf*ones(3*I+J,1); G; NaN*ones(J,1); Inf*ones(K_Flood,1)];        
    
    lb = full(repmat(lb,T_opt,1));
    ub = full(repmat(ub,T_opt,1));

    for j=1:J
        ub((0:T_opt-1)*K+3*J+3*I+j) = Irr(Month_opt,j);
        % Fix irrigation variables, comment out for stage 1.
        lb((0:T_opt-1)*K+3*J+3*I+j) = Irr(Month_opt,j);
    end

%%  Constraints bounds:cl < [M*x ; N*x] < cu
    
    w = (CatchmA'./Area)*Qm_opt';
    w = w(:);
    w = -kron(eye(T_opt),D)*w; 
    
    % System without irrigation inequality constraints, stage 2:
    if K_Flood > 0  % If floods are included in the optimization
        maxStor = lev_Flood./Conv(seg_Flood);    
        d = repmat([zeros(J,1) ; maxStor'],T_opt,1);   
        if ~isnan(consIncluded)
            d = repmat([zeros(length(consIncluded),1) ; maxStor'],T_opt,1);
        end
    else
        d = repmat(zeros(J,1),T_opt,1);  
        if  ~isnan(consIncluded)
            d = repmat(zeros(length(consIncluded),1),T_opt,1);
        end
    end
    
% % System with irrigation inequality constraints, stage 1:
% if K_Flood > 0  % If floods are included in the optimization
%     maxStor = lev_Flood./Conv(seg_Flood);    
%     d = repmat([zeros(2*J,1) ; maxStor'],T_opt,1);  
%     if ~isnan(consIncluded)
%       d = repmat([zeros(2*length(consIncluded),1) ; maxStor'],T_opt,1);
%     end
% else
%     d = repmat(zeros(2*J,1),T_opt,1);  
%     if  ~isnan(consIncluded)
%       d = repmat(zeros(length(consIncluded),1),T_opt,1);
%     end
% end
    
    cl = [w ; -Inf*ones(length(d),1)];
    cu = [w ; d];          

end
    
    



