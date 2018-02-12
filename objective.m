%   Objective function

function [ obj ] = objective(x,auxdata)

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
    G  = auxdata.ResDetails.G; 
    Irr = mean(auxdata.ResDetails.Irr);    
    
    % System details
    J = auxdata.system.J;
    I = auxdata.system.I;
    Q2 = auxdata.system.Q2;
    K = auxdata.system.K;
    
%%  Hydropower and irrigation term in the objective    
    
    % Include the hydropower and irrigation only for considered reservoir
    if ~isnan(consIncluded)
        % Hydropower term
        H1 = sparse(2*J+3*I+consIncluded,consIncluded,...
            1./(T_opt.*G(consIncluded).*S(consIncluded)),K,K);
        H = kron(speye(T_opt),H1); 

        % Irrigation term
        f1 = sparse(1,(3*J+3*I+consIncluded),...
            1./(T_opt.*Irr(consIncluded)),1,K);    
        f = repmat(f1,1,T_opt);
    else
        % Hydropower term
        H1 = sparse(2*J+3*I+(1:J), 1:J ,1./(T_opt.*G.*S),K,K);
        H = kron(speye(T_opt),H1); 

        % Irrigation term
        f1 = sparse(ones(J,1),(3*J+3*I+(1:J)),1./(T_opt.*Irr),1,K);    
        f = repmat(f1,1,T_opt);
    end
    
%%  Flood term in the objective function
%   This term is zero when flood are not included in the objective, when 
%   K_Flood = 0.
    
    Pen_ConsFlood = sparse(ones(K_Flood,1),4*J+3*I+(1:K_Flood),flNorm,1,K);    
    Pen_ConsFlood = repmat(Pen_ConsFlood,1,T_opt);
    
%%  Spillway term in the objective function
%   This term is zero when no spillway term is included in the
%   objective, when constrSpill = 0.

    Pen_ConsSpill = 0;  
 
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
        g = Mnl*x -  Qin + Smax;                        

        % Pen_ConsSpill := 
        %   Spillway*(Smax - Storage + Spillway + Conduit - Inflow)
        bnl = sparse(1:J,J+3*I+(1:J),ones(J,1),J,K);    
        bnl = kron(eye(T_opt),bnl);
        Pen_ConsSpill = (bnl*x).*g;                   
      
        % Include the spillway constraint only for considered reservoir 
        if  ~isnan(consIncluded)
            Pen_ConsSpill = ...
                Pen_ConsSpill([consIncluded+(0:(T_opt-1))*J],:);      
        end    
        
        Pen_ConsSpill = sum(Pen_ConsSpill);
                
    end     
    
%%  Objective function
%   beta = [weight HP, weight Irr], 
%   rho = penalty for the spillway
%   rho2 =   penalty for floods
%   Hydropower = x'*H*x, irrigation = f*x

    obj = (-1*(beta(1)*x'*H*x + beta(2)*f*x ...
        - rho*Pen_ConsSpill - rho2*Pen_ConsFlood*x));
       
end