%%  Optimization and river simulation for stage 1

    clear all
    close all hidden
    clc 
    
    cd(fileparts(which(mfilename)));
        
    DateTime = clock;

%%  Settings 

%   Optimization parameters
    % Start timestep of optimization
    T_start = 1;                    
    % Length of the optimization horizon (days)
    T_opt = 365;         
    % Length of the prediction horizon (days)
    T_pred = 3;               
    % Number of timesteps the horizon is shifted each time
    T_shift = 1;  
    % Number of times the horizon is shifted 
    T_mov = 1500;          
    % Length of timestep (s)    
    dt = 24*3600;         
    % Tolerance level for the overall NLP error IPOPT
    T_tol = 100;               
    % Objective function scaling factor
    ObjScaling = T_opt*10^5; 
    % Spillway constraint, 0 = not included, 2 = included
    constrSpill = 0;                
    % If constrSpill == 2, how to the spillway constraint
    % 0 = exact, 2 = penalty method
    constrMethod = 0;   
    % Penalty parameter  
    rho =  0;     
    
    % Flood constraint
    % Include flood constraint (1 = yes, 0 = no)
    constrFlood = 0;   
    % Number segments with flood constraint
    K_Flood = 0;   
    % Segments with flood constraint
    seg_Flood = 7;             
    % Water level resulting in floods
    lev_Flood = 12;  
    % Penalty for flood constraint   
    rho2 = 0;              
    % Normalization factor for flood term
    flNorm = 1/(T_opt*10^4*37);
     
    % Make structure for settings of flood constraint
    flood.K_Flood = K_Flood;
    flood.seg_Flood = seg_Flood;
    flood.lev_Flood = lev_Flood;
    flood.rho2 = rho2;
    flood.flNorm = flNorm;
        
    % Weight vector objective function
    % beta(1) = weight HP, beta(2) = weight irrigation
    beta = [1 0];     
    % Pick which reservoirs to include in the model
    Index = (1:21);                 
    % Number of timesteps for which results are saved
    T = T_mov;      
    
    optimization.T_start = T_start;
    optimization.T_opt = T_opt;
    optimization.T_pred = T_pred;
    optimization.T_shift = T_shift; 
    optimization.T_mov = T_mov;
    optimization.T_tol = T_tol;
    optimization.T = T;
    optimization.dt = dt;
    optimization.constrSpill = constrSpill;
    optimization.constrMethod = constrMethod;
    optimization.constrFlood = constrFlood; 
    optimization.rho = rho;
    optimization.beta = beta;
    optimization.Index = Index;
    optimization.flood = flood;  
    optimization.varboundsFixed = [];
    optimization.xFixed = [];
    optimization.consIncluded = [];
    optimization.ObjScaling = ObjScaling;
    
%%  Run MPC stage 1
     
    optimization.MPC = 1;
    [system,ResDetails,RivDetails] = systemDynamics(optimization);  
    
    %   Save results and scripts -----------------------------------------
    
        Run = '1';
        saveTo = ['results\run' , Run , '\'];
        
        % Save results     
        save([saveTo, 'run' , Run , '_' , num2str(DateTime(3)) ,'-' , ...
            num2str(DateTime(2)), '-', num2str(DateTime(4)) ,'h', ...
            num2str(DateTime(5)) ,'.mat' ],'optimization','system', ...
        'ResDetails','RivDetails');     

      % Save current file
        FileNameAndLocation=[mfilename(['fullpath'])];
        BackupNameAndLocation=[[saveTo , ['runStage1_' ...
            num2str(DateTime(3)) '-' num2str(DateTime(2)) '-'...
            num2str(DateTime(4)) 'h' num2str(DateTime(5))  ]]];
        newbackup=sprintf('%s_run%s.m',BackupNameAndLocation,Run);
        currentfile=strcat(FileNameAndLocation, '.m');
        copyfile(currentfile,newbackup);

        % Save runOptimization.m   
        FileNameAndLocation=[cd , '\runOptimization_Stage1'];
        BackupNameAndLocation=[saveTo , ['runOptimization_Stage1_' ... 
            num2str(DateTime(3)) '-' num2str(DateTime(2)) '-'...
            num2str(DateTime(4)) 'h' num2str(DateTime(5)) ]];
        newbackup=sprintf('%s_run%s.m',BackupNameAndLocation,Run);
        currentfile=strcat(FileNameAndLocation, '.m');
        copyfile(currentfile,newbackup);   
    %   ------------------------------------------------------------------     

    J = system.J;
    K_Riv = system.K_Riv;       
    A = system.A;
    B = system.B;
    C = system.C;
    D = system.D;    
    
    CatchmA = ResDetails.CatchmA;
       
    Status = NaN(J,T_mov);
    Iter = NaN(J,T_mov);
    CPU = NaN(J,T_mov);
    Xstor = NaN(J,T);
    Xspill = NaN(J,T);
    U = NaN(2*J,T);
    UB1 = NaN((4*J+K_Riv),T);

for j = 1:21  
    
    % Waitbar
    h = waitbar(0,['Please wait run ' num2str(Run) '  (' ...
        num2str(j)  '/' num2str(21) ') ...']);
              
        optimization.Index = Index(j); 
        optimization.MPC = 2;
        [system,ResDetails,RivDetails] = systemDynamics(optimization);
    
        [results1,forcing] = runOptimization_Stage1(optimization,...
            system,ResDetails,RivDetails);

        % For saving results stage 1 for this reservoir
        X2 = results1.X;      
        UB2 = results1.UB;
        Status(j,:) = results1.Status;
        Iter(j,:) = results1.Iter;
        CPU(j,:) = results1.CPU ;
        
        % For running the complete system        
        Xstor(j,:) = X2(1+4*((1:T)-1));            
        Xspill(j,:) = X2(2+4*((1:T)-1));             
        U(j,:) = X2(3+4*((1:T)-1));                 
        U(J+j,:) = X2(4+4*((1:T)-1));               

        Xstor(j,:) =  X2(1+4*((1:T)-1));           
        Xspill(j,:) = X2(2+4*((1:T)-1));            
        U(j,:) = X2(3+4*((1:T)-1));                
        U(J+j,:) =  X2(4+4*((1:T)-1));                
        
        UB1(j,:) = UB2(1+4*((1:T)-1));
        UB1(J+K_Riv+j,:) = UB2(2+4*((1:T)-1));
        UB1(2*J+K_Riv+j,:) = UB2(3+4*((1:T)-1));
        UB1(3*J+K_Riv+j,:) = UB2(4+4*((1:T)-1));
                
    close all hidden
    
end    
    
    Qm = forcing.Qm;   
    Area = forcing.Area;
    W = (CatchmA'./Area)*Qm((1:T))';  
       
    x_init = zeros(K_Riv,1);
    x = x_init;
   
    Xriv = zeros(K_Riv,T);
    Xriv(:,1) = x;
    
% Simulate river for operations of entire system    
for k = 1:(T-1)
    
    Xriv(:,k+1) = B(J+(1:K_Riv),J+(1:K_Riv))\...
        (A(J+(1:K_Riv),J+(1:K_Riv))*Xriv(:,k) ...
        + A(J+(1:K_Riv),J+K_Riv+(1:J))*Xspill(:,k) ...
        + C(J+(1:K_Riv),:)*U(:,k) + D(J+(1:K_Riv),:)*W(:,k));

end
    
    X1 = [ Xstor ; Xriv ; Xspill ; U ];
    X1 = X1(:);
    UB1 = UB1(:);

   results.X = X1;
   results.UB = UB1;
   results.Status = Status;    
   results.Iter = Iter;
   results.CPU = CPU;   
   
   % Save results stage 1 for entire system
   save([saveTo, 'run' , Run , '_' , num2str(DateTime(3)) ,'-' , ...
        num2str(DateTime(2)), '-', num2str(DateTime(4)) ,'h', ...
        num2str(DateTime(5)) ,'.mat' ],'results','forcing','-append'); 
 
 