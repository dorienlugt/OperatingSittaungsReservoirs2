%%  Stage 2 optimization

    clear all
    close all hidden
    clc
    
    cd(fileparts(which(mfilename)));
    
    DateTime = clock;
    
%%  Settings, results and parameters from stage 1   
        
%   Save results and scripts -----------------------------------------

        Run = '2';   

        % Load stage 1 results
        output1 = load(['results\run53\run53.mat']);
        saveTo = ['results\run' , Run , '\'];
          
        % Save current file
        FileNameAndLocation=[mfilename(['fullpath'])];
        BackupNameAndLocation=[[saveTo , ['runStage2_'...
            num2str(DateTime(3)) '-' num2str(DateTime(2)) '-'...
            num2str(DateTime(4)) 'h' num2str(DateTime(5))  ]]];
        newbackup=sprintf('%s_run%s.m',BackupNameAndLocation,Run);
        currentfile=strcat(FileNameAndLocation, '.m');
        copyfile(currentfile,newbackup);

        % Save runOptimization.m   
        FileNameAndLocation=[cd , '\runOptimization_Stage2'];
        BackupNameAndLocation=[saveTo , ['runOptimization_Stage2_'  ...
            num2str(DateTime(3)) '-' num2str(DateTime(2)) '-'...
            num2str(DateTime(4)) 'h' num2str(DateTime(5)) ]];
        newbackup=sprintf('%s_run%s.m',BackupNameAndLocation,Run);
        currentfile=strcat(FileNameAndLocation, '.m');
        copyfile(currentfile,newbackup);   
    
%   ------------------------------------------------------------------     

    % Load settings from stage 1
    optimization = output1.optimization;
    system = output1.system;
    ResDetails = output1.ResDetails;
    RivDetails = output1.RivDetails;
    forcing = output1.forcing;
    results = output1.results; 
    
    % Optimization horizon (days) 
    T_opt = 30;                                   
    optimization.T_opt = T_opt;
    % Number of times horizon is shifted
    T_mov = 10; %365;             
    % Number of timesteps horizon is shifted each shift
    T_shift = optimization.T_shift;        
    % Start on April 1st 2011
    T_start = 455;       
    % Objective function scaling factor
    ObjScaling = T_opt*10^5; 
    
    % Number of reservoirs      
    J = system.J;     
    % Number of river sections     
    I = system.I;   
    % Number of variables in stage 1
    K1 = system.K;                                  
    % Conversion factors for storage to water level in river sections
    Conv = RivDetails.Conv;                        
    
    % The reservoirs included in the system  
    Index = optimization.Index;   
    % For each catchment the river section to which it is connected
    ResCatchm = ResDetails.ResCatchm;
                  
    % Results of stage 1
    X1 = results.X;                                
    X1 = reshape(X1,K1,length(X1)/K1);  
    
%%   Settings and parameters in stage 2    
        
    % Flood constraint
    % Include flood constraint (1 = yes, 0 = no)
    constrFlood = 1;   
    % Number segments with flood constraint
    K_Flood = 1;   
    % Segments with flood constraint
    seg_Flood = 7;             
    % Water level resulting in floods
    lev_Flood = 12;  
    % Penalty for flood constraint   
    rho2 = 0;              
    % Normalization factor for flood term
    flNorm = 1/(T_opt*10^4*37);
    
    % Make structure for settings of flood constraint
    optimization.constrFlood = constrFlood;      
    optimization.flood.K_Flood = K_Flood;          
    optimization.flood.seg_Flood = seg_Flood;       
    optimization.flood.lev_Flood = lev_Flood;      
    optimization.flood.rho2 = rho2;           
    optimization.flood.flNorm = flNorm;
    
    % Initial state for stage 2 from results of stage 1
    if  K_Flood == 1
        xFixed = [ X1;...
            max(X1(J+3*seg_Flood,:) - 0.9*lev_Flood/Conv(seg_Flood),0) ];         
    else
        xFixed = X1; 
    end
    optimization.xFixed = xFixed(:); 
    optimization.T_start = T_start;  
    optimization.T_mov = T_mov;
    optimization.ObjScaling = ObjScaling;
        
    % Spillway constraint, 0 = not included, 2 = included
    optimization.constrSpill = 2;              
    % If constrSpill == 2, how to the spillway constraint
    % 0 = exact, 2 = penalty method
    optimization.constrMethod = 2;          
    optimization.MPC = 1;
    % Weight vector objective function
    % beta(1) = weight HP, beta(2) = weight irrigation
    optimization.beta = 1*[1 0];   
    
    optimization.rhoIterates = [0 10^(-15) 10^(-14)]; 
    LrhoIt =  length(optimization.rhoIterates);
    
%   Stage 2 optimization
    % All reservoirs upstream of the river section with flood constraint
    for j =  Index(ResCatchm<=seg_Flood)            
 
        % Save results in file:    
        save([saveTo, 'run' num2str(Run) '_Res' num2str(j)  ...
            '_' num2str(DateTime(3)) '-' num2str(DateTime(2)) '-' ...
            num2str(DateTime(4)) 'h' num2str(DateTime(5)) ],'output1'); 
 
        Status = NaN(LrhoIt,T_mov);
        Iter = NaN(LrhoIt,T_mov);
        CPU = NaN(LrhoIt,T_mov);     
        X2 = NaN(LrhoIt,((T_mov-1)*T_shift+T_opt)*(4*J+3*I+K_Flood)); 
                
        % Include constraints for the reservoir adjusted in this iteration
        optimization.consIncluded = j;                  

        [system,ResDetails,RivDetails] = systemDynamics(optimization);  
        % Number of variables for stage 2
        K2 = system.K;
        K_Riv = system.K_Riv;
      
        % The variables fixed for this iteration of stage 2
        varboundsFixed = ((1:T_opt)-1)'*K2 + [Index(Index~=j) ...
            J+K_Riv+Index(Index~=j) 2*J+K_Riv+Index(Index~=j) ...
            3*J+K_Riv+Index(Index~=j)];  
        varboundsFixed = varboundsFixed';
        optimization.varboundsFixed = varboundsFixed(:);
                      
        % Iterate over all values of the penalty method 
        % for the spillway term  
        for m = 1:LrhoIt

            h = waitbar(0,['Please wait Res '...
               num2str(j) ' (' num2str(m)  '/' num2str(LrhoIt) ') ...']);

            optimization.rho = optimization.rhoIterates(m);         
            [results,forcing] = runOptimization_Stage2(optimization, ...
                system,ResDetails,RivDetails);

            X2(m,1:length(results.X)) = results.X';
            optimization.xFixed(((T_start-1)*K2) ...
                +(1:length(X2(m,:)))) = X2(m,:);
            optimization.xFixed((((T_start-1)*K2) ...
                +length(X2(m,:))+1):end) = NaN;

            Status(m,:) = results.Status;
            Iter(m,:) = results.Iter;
            CPU(m,:) = results.CPU; 

            close all hidden

        end
    
        results2.X = X2;
        results2.UB = results.UB;
        results2.Status = Status;
        results2.Iter = Iter;
        results2.CPU = CPU;   
    
     	optimization2 = optimization;
        forcing2 = forcing;
        ResDetails2 = ResDetails;
        RivDetails2 = RivDetails;
        system2 = system;
        
        save(['results\run' num2str(Run) '\run' num2str(Run) ...
            '_Res' num2str(j)  '_' num2str(DateTime(3)) '-'...
            num2str(DateTime(2)) '-' num2str(DateTime(4)) 'h'...
            num2str(DateTime(5)) '.mat'],'optimization2','ResDetails2',...
            'RivDetails2','system2','results2','forcing2','-append');
        
    end          
       
    
    

 