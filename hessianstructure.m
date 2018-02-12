%   Structure of the Hessian matrix of the Lagrangian 

function hstr = hessianstructure(auxdata)      

    % Settings for the optimization 
    T_opt = auxdata.optimization.T_opt;

    % System details
    J = auxdata.system.J;
    I = auxdata.system.I; 
    K = auxdata.system.K;      
       
    hstr2 = sparse([1:J 1:J 1:J J+3*I+(1:J) J+3*I+(1:J) 2*J+3*I+(1:J)],...
        [1:J J+3*I+(1:J) 2*J+3*I+(1:J) J+3*I+(1:J) 2*J+3*I+(1:J)...
        2*J+3*I+(1:J)],1,K,K);
    hstr = kron(eye(T_opt),hstr2); 

    hstr = sparse(hstr);
    hstr = hstr'; 
    
end