%   Gradient of the objective function

function grad = gradient(x,auxdata)

    % Setting for the optimization
    T_opt = auxdata.optimization.T_opt;
    beta = auxdata.optimization.beta;             
    constrSpill = auxdata.optimization.constrSpill;
    constrMethod = auxdata.optimization.constrMethod; 
    consIncluded = auxdata.optimization.consIncluded;     
    K_Flood = auxdata.optimization.flood.K_Flood;
    rho = auxdata.optimization.rho;
    rho2 = auxdata.optimization.flood.rho2;
    flNorm = auxdata.optimization.flood.flNorm;
    
    % Forcing data
    Area = auxdata.forcing.Area;
    Qm_opt = auxdata.forcing.Qm_opt;
    
    % Reservoir details 
    CatchmA = auxdata.ResDetails.CatchmA;
    S  = auxdata.ResDetails.S; 
    Gcond  = auxdata.ResDetails.G; 
    Irr = mean(auxdata.ResDetails.Irr);

	% System details
	J = auxdata.system.J;
    I = auxdata.system.I;  
    K = auxdata.system.K;  
    Q2 = auxdata.system.Q2;
    
%%  Gradient of the hydropower and irrigation term in the objective    
    
    % Include the hydropower and irrigation only for considered reservoir
    if ~isnan(consIncluded)
        % Gradient of the hydropower term
        G1 = sparse(consIncluded, 2*J+3*I+consIncluded,...
            1./(T_opt*Gcond(consIncluded).*S(consIncluded)),K,K);
        G = kron(speye(T_opt),G1) ;

        G2 = sparse(2*J+3*I+consIncluded,consIncluded,...
            1./(T_opt*Gcond(consIncluded).*S(consIncluded)),K,K);
        G = G + kron(speye(T_opt),G2);        
        
        % Gradient of the irrigation term
        g1 = sparse(1,(3*J+3*I+consIncluded),...
            1./(T_opt*Irr(consIncluded)),1,K);    
        g = repmat(g1,1,T_opt);
    else
        % Gradient of the hydropower term
        G1 = sparse(1:J, 2*J+3*I+(1:J),1./(T_opt*Gcond.*S),K,K);
        G = kron(speye(T_opt),G1) ;

        G2 = sparse(2*J+3*I+(1:J),1:J,1./(T_opt*Gcond.*S),K,K);
        G = G + kron(speye(T_opt),G2);        

        % Gradient of the irrigation term
        g1 = sparse(ones(J,1),(3*J+3*I+(1:J)),1./(T_opt*Irr),1,K);    
        g = repmat(g1,1,T_opt);
    end
    
%%  Gradient of the flood term in the objective function
%   This term is zero when flood are not included in the objective, when 
%   K_Flood = 0.
    
    gradPen_ConsFlood = sparse(ones(K_Flood,1),4*J+3*I+(1:K_Flood),...
        flNorm,1,K);    
    gradPen_ConsFlood = repmat(gradPen_ConsFlood,1,T_opt);
    
%%  Gradient of the spillway term in the objective function
%   This term is zero when no spillway term is included in the
%   objective, when constrSpill = 0.
    
    gradPen_ConSpill = 0;         
	  
    if  constrSpill == 2 && constrMethod == 2     
        
        % Inflows (Qin) to the reservoir(s)
        w = (CatchmA'./Area)*Qm_opt';     
        Qin = Q2*w;
        Qin = Qin(:);
        % Storage capacities (Smax) for the reservoir(s) 
        Smax = repmat(S,T_opt,1);

        % g := 
        %   Smax - Storage + Conduit + Spillways - Inflow
        Mnl = sparse([1:J 1:J 1:J],[1:J J+3*I+(1:J) 2*J+3*I+(1:J)],...
            [-ones(J,1) ones(J,1) ones(J,1)],J,K); 
        Mnl = kron(eye(T_opt),Mnl);

        % gSpill := 
        %   Smax - Storage + Conduit + Spillways - Inflow 
        gSpill = Mnl*x -  Qin + Smax;     
        % bnl :=
        %   Spillway
        bnl = sparse(1:J,J+3*I+(1:J),ones(J,1),J,K);  
        bnl = kron(eye(T_opt),bnl)*x;                        

        gradPen_ConSpill = sparse(0,0);

        for t = 1:T_opt

            Mat1 = sparse((1:J),(1:J),bnl(J*(t-1)+(1:J)),J,J);
            Mat2 = sparse((1:J),(1:J),...
                    gSpill(J*(t-1)+(1:J))+bnl(J*(t-1)+(1:J)),J,J);
            Mat = [-Mat1 sparse(J,3*I) Mat2 Mat1 sparse(J,J) ...
                    sparse(J,K_Flood)];  
            gradPen_ConSpill = blkdiag(gradPen_ConSpill,Mat);

        end    
        
        % Include the spillway constraint only for considered reservoir 
        if  ~isnan(consIncluded)

            gradPen_ConSpill = ...
                gradPen_ConSpill(consIncluded+(0:T_opt-1)*J,:);

        end

        gradPen_ConSpill = sum(gradPen_ConSpill);

    end    
    
%%  Gradient vector of the objective function
%   beta = [weight HP, weight Irr], 
%   rho = penalty for the spillway
%   rho2 = penalty for floods
%   Hydropower = (G*x)', irrigation = g

    grad = -1*(beta(1)*(G*x)' + beta(2)*g ...
        - rho*gradPen_ConSpill - rho2*gradPen_ConsFlood);    
   
end