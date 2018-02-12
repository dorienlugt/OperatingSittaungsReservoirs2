%  Initial solution for optimization at the first timestep in stage 1

function x_init = initial(auxdata)

    % All reservoirs are filled up to 64% of the storage capacity, similar
    % to reservoir 3 on March 1st, 2010. In initial solution, no flow is 
    % released through the spillway, the conduits release water at 50% of 
    % the conduit capacity, irrigation offtake is equal to the demand. 
    
    T_opt = auxdata.optimization.T_opt;
    
    S = auxdata.ResDetails.S;      
    G = auxdata.ResDetails.G;   
        
    J = auxdata.system.J;
    I = auxdata.system.I;
    K = auxdata.system.K;
    
    [~,ub,~,~] = bounds(auxdata);       

    x_init = repmat([(231088824/360000000)*S; ones(3*I,1); ...
        zeros(J,1) ; 0.5*G ; NaN(J,1) ],T_opt,1); 
    
    for j=1:J
        x_init((0:T_opt-1)*K+3*J+3*I+j) = ub((0:T_opt-1)*K+3*J+3*I+j);
    end
    
end


