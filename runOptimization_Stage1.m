
function [results,forcing] = ...
    runOptimization(optimization,system,ResDetails,RivDetails)

    T_start = optimization.T_start;
    T_opt = optimization.T_opt;
    T_pred = optimization.T_pred;
    T_mov = optimization.T_mov;
    T_shift = optimization.T_shift;
    T_tol = optimization.T_tol;
    xFixed = optimization.xFixed;
    varboundsFixed = optimization.varboundsFixed;
    optimization.flood.K_Flood;     
    ObjScaling = optimization.ObjScaling;
    
    M = system.M;
    K = system.K;
    
    [Qm,Qseas,Day,Month,Year,Area] = ...
        forcingdata_Seasonality(optimization,ResDetails);
    Month_opt = Month(1:T_opt);
    Qm_opt = Qseas(1:T_opt);        
    Qm_opt(1:T_pred) = Qm(1:T_pred); 
    
    opts.auxdata.forcing.Day = Day;
    opts.auxdata.forcing.Month = Month;
    opts.auxdata.forcing.Year = Year;
    opts.auxdata.forcing.Qm = Qm;
    opts.auxdata.forcing.Qseas = Qseas;
    
    opts.auxdata.forcing.Month_opt = Month_opt;
    opts.auxdata.forcing.Qm_opt = Qm_opt;
    opts.auxdata.forcing.Area = Area;
    opts.auxdata.optimization = optimization;
    opts.auxdata.system = system;
    opts.auxdata.ResDetails = ResDetails;
    opts.auxdata.RivDetails = RivDetails;

    funcs.objective = @objective;                  
    funcs.gradient = @gradient;                    
    funcs.constraints = @constraints;               
    funcs.jacobian = @jacobian;                     
    funcs.jacobianstructure = @jacobianstructure;   
    funcs.hessian = @hessian;                       
    funcs.hessianstructure = @hessianstructure;       

    if isempty(xFixed)    
        x_init = initial(opts.auxdata);        
    else
        x_init = xFixed(K*(T_start-1)+(1:T_opt*K));
    end
    
    [~,~,cl,~] = bounds(opts.auxdata);

    % Settings for IPOPT
    opts.ipopt.print_level  = 1;  
    opts.ipopt.max_iter = 3000;
    opts.ipopt.tol = T_tol;
    opts.ipopt.constr_viol_tol = 100;
    opts.ipopt.compl_inf_tol = 100;
    opts.ipopt.dual_inf_tol = 1;
    opts.ipopt.acceptable_constr_viol_tol = 10^3;
    opts.ipopt.acceptable_compl_inf_tol = 10^3;
    opts.ipopt.obj_scaling_factor = ObjScaling;

    cl_extra = zeros(length(cl),1);  
    [Mr,~] = size(M);
    cl_extra(1:Mr/T_opt) = -M(Mr/T_opt+(1:Mr/T_opt),...
        1:(length(x_init)/T_opt))*x_init(1:(length(x_init)/T_opt));
    cu_extra = cl_extra;
    
    X =zeros(0); 
    UB = zeros(0);
    
    Status = NaN(T_mov,1);
    Iter = NaN(T_mov,1);
    CPU = NaN(T_mov,1); 
    
for k = 1:T_mov
    
        Qm_opt = Qm((k-1)*T_shift+(1:T_opt));       
        Qm_opt(1:T_pred) = Qseas((k-1)*T_shift+(1:T_pred)); 
    
    
        Month_opt = Month((k-1)*T_shift+(1:T_opt));
        opts.auxdata.forcing.Month_opt = Month_opt;
        opts.auxdata.forcing.Qm_opt = Qm_opt;
           
        [lb,ub,cl,cu] = bounds(opts.auxdata);      
        
        cl = cl + cl_extra;
        cu = cu + cu_extra;    
        
        lb(varboundsFixed) = xFixed((T_start-1+(k-1)*T_shift)*K + ...
            varboundsFixed);
        ub(varboundsFixed) = xFixed((T_start-1+(k-1)*T_shift)*K + ...
            varboundsFixed);
        
        opts.cl = cl;           % constraints lower bounds
        opts.cu = cu;           % constraints upper bounds
        opts.lb = lb;           % variable lower bounds
        opts.ub = ub;           % variable upper bounds  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
        [x,info] = ipopt(x_init,funcs,opts);    
        
        Status(k) = info.status;
        Iter(k) = info.iter;
        CPU(k) = info.cpu;

        if  isempty(xFixed)    
       
            x_init = [x((T_shift*K+1):end) ; ...
                x((T_opt-T_shift)*K + (1:T_shift*K))];  
     
        else
            x_init = [x((T_shift*K+1):end) ;  ...
                xFixed(((T_start-1)+(k-1)*T_shift)*K + ...
                (T_opt-T_shift)*K + (1:T_shift*K)) ]; 
        end

        cl_extra = [ -M(Mr/T_opt+(1:Mr/T_opt),...
            1:(K))*x((T_shift-1)*K+(1:K)); zeros(length(cl)-Mr/T_opt,1)] ;   
        cu_extra = cl_extra;   
                
        UB = [UB ; ub(1:T_shift*length(ub)/T_opt)];                       
        X = [X ; x(1:T_shift*length(x)/T_opt)];
           
end

        UB = [UB ; ub((T_shift*length(ub)/T_opt+1):end)];                
        X = [X ; x((T_shift*length(x)/T_opt+1):end)];
            
        results.X = X;
        results.UB = UB;
        results.Status = Status;
        results.Iter = Iter;
        results.CPU = CPU;
        
        forcing = opts.auxdata.forcing;

end

