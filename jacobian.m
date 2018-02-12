%   Jacobian matrix of the constraints

function jac = jacobian( ~ ,auxdata)        
    
    N = auxdata.system.N;   % Linear inequality constraints
    M = auxdata.system.M;   % Linear equality constraints
    
    jac = [M ; N ];  
   
 
end