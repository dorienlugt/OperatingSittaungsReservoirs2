%   Constraints

function    cons = constraints(x,auxdata)
    
    N = auxdata.system.N;   % Linear inequality constraints
    M = auxdata.system.M;   % Linear equality constraints

    cons = full([M*x ; N*x ]);

end
