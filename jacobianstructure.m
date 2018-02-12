%   Structure of the Jacobian matrix of the constraints

function jacstr = jacobianstructure(auxdata)
    
    N = auxdata.system.N;   % Linear inequality constraints
    M = auxdata.system.M;   % Linear equality constraints 
    
    jacstr = spones([M ; N]);     

    
end     



