%   Hessian of the Lagrangian

function Hessian = hessian(~,sigma,~,auxdata)

    % Setting for the optimization
    T_opt = auxdata.optimization.T_opt;
    beta = auxdata.optimization.beta;            
    constrSpill = auxdata.optimization.constrSpill;    
    constrMethod = auxdata.optimization.constrMethod;
    consIncluded = auxdata.optimization.consIncluded; 
    rho = auxdata.optimization.rho;
    
    % Reservoir details 
    S = auxdata.ResDetails.S;    
    G = auxdata.ResDetails.G;   
    
    % System details
    J = auxdata.system.J;
    I = auxdata.system.I;  
    K = auxdata.system.K;  
    
%%  Second order derivative of the hydropower term in the objective    
      
    % Include the hydropower term only for considered reservoir
    if  ~isnan(consIncluded)
        Hobj = sparse(consIncluded,2*J+3*I+consIncluded,...
            1./(T_opt*G(consIncluded).*S(consIncluded)),K,K);
    else
        Hobj = sparse(1:J,2*J+3*I+(1:J),1./(T_opt*G.*S),K,K);
    end
        Hobj =  -beta(1)*kron(eye(T_opt),Hobj); 

%%  Second order derivative of the spillway term in the objective
%   This term is zero if spillway term is not included in the objective,
%   when constrSpill = 0.

    HobjPen = sparse(T_opt*K,T_opt*K);

        if  constrSpill == 2 && constrMethod == 2                
            
            % Include the hydropower term only for considered reservoir     
            if  ~isnan(consIncluded)
                j = consIncluded;
                HobjPen = sparse([j J+3*I+j J+3*I+j],... 
                    [J+3*I+j J+3*I+j 2*J+3*I+j],[-1 2 1],K,K);
                HobjPen = -rho*kron(eye(T_opt),HobjPen);   

            else                       

                HobjPen = sparse([1:J J+3*I+(1:J) J+3*I+(1:J)],...
                    [J+3*I+(1:J) J+3*I+(1:J) 2*J+3*I+(1:J)],...
                    [-ones(J,1) 2*ones(J,1) ones(J,1)],K,K);
                HobjPen = -rho*kron(eye(T_opt),HobjPen); 

            end

        end
    
%%  Hessian matrix of the Lagrangian

    Hessian = sigma*(Hobj + HobjPen);   
    Hessian = Hessian';       
    
end
