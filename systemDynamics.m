%%  Model of the system dynamics

function [system,ResDetails,RivDetails] = systemDynamics(optimization)

    T_opt = optimization.T_opt;
    Index = optimization.Index;  
    MPC = optimization.MPC;
    K_Flood = optimization.flood.K_Flood;           
    seg_Flood = optimization.flood.seg_Flood;      
    consIncluded = optimization.consIncluded;      
    
    [ResDetails] = reservoirDetails(optimization);
    [RivDetails] = riverDetails(optimization);
        
    ResCatchmA = ResDetails.ResCatchmA;
    ResCatchm = ResDetails.ResCatchm;
    CatchmA = ResDetails.CatchmA;
    
    c1 = RivDetails.c1;
    c2 = RivDetails.c2;
    c3 = RivDetails.c3;            
            
    if  MPC == 1                % MPC stage 1   
        
        % Number of reservoirs 
        J = length(Index);     
        % Number of river segments 
        I = 17;           
        
        % Matrix JxI with 1 at (j,i) if reservoir j is in catchment i
        Q1 = sparse(Index,ResCatchm,1);            
        % Matrix JxI with at (j,i) the reservoir catchment area j 
        Q3 = sparse(Index,ResCatchm,ResCatchmA); 
        % Ratio of reservoir area j of catchment area i 
        Q2  = Q3./CatchmA;                           
        
    else
        if  MPC == 2            % MPC stage 1            
            
            % Number of reservoirs 
            J = 1;
            % Number of river segments 
            I = 0;
            
            % Ratio of reservoir area j of catchment area i        
            Q2 = ResCatchmA/CatchmA;  
            
        end
    end
    
    % Number of state veriables for river segments
    K_Riv = 3*I;
    % Number of state variables
    K = 4*J + K_Riv + K_Flood;
            
%%  Create matrices A, B, C, D 
%   with B x^(k+1) = A x^(k) + C u^(k) + D w^(k), where     
%   x^(k) is the vector of state variables at timestep k
%   u^(k) is the vector of control variables at timestep k
%   w^(k) is the vector of external disturbances at timestep k

    %   Equality constraints for reservoir dynamics  
    %   Storage in the reservoir minus flow through the conduit
    B1 = speye(J); 
    A1 = speye(J); 
    C1  = [-speye(J) sparse([],[],[],J,J)];
    %   Inflow to the reservoir from the reservoir catchment
    D1 = Q2;

    A = sparse([],[],[],0,0);
    B = sparse([],[],[],0,0);
    C = sparse([],[],[],0,0);
    D = sparse([],[],[],0,0);  
    A4 = sparse([],[],[],0,0);                       

    % Equality constraints for river segments
    for  i = 1:I
        
        %  Muskingum routing equations
        A2 = sparse([2 2 3 3 3],[1 2 1 2 3],...
            [c1(i) c3(i) 1/2 -1/2 1],3,3);
        B2 = sparse([1:3 2 3 3],[1:3 1 1 2],...
            [ones(1,3) -c2(i) -1/2 1/2],3,3);

        % Flow through the conduit - irrigation offtake to river section     
        C2 = [Q1(:,i)' -Q1(:,i)' ; sparse([],[],[],2,2*J)];    
        % Water from uncontrolled part of the catchment to the river   
        D2 = sparse(1,i,(1-sum(Q2(:,i))),3,I);   
        
        C = [C ; C2];
        D = [D ; D2];
        A = blkdiag(A,A2);
        B = blkdiag(B,B2);

        % Water from the spillway to the river section
        A3 = sparse(ones(J,1),1:J,Q1(:,i),3,J);                 
        A4 = [A4; A3];
        
    end   
    
   % Add outflow from previous section to new section
   B = B + sparse(3*(2:I)-2,3*(1:(I-1))-1,-ones(I-1,1),3*I,3*I); 
   
   A = blkdiag(A1,A);
   B = blkdiag(B1,B);        

   % Storage minus outflow through the spillway
   A = [A , [-speye(J) ; A4]];
   
   B = [B , sparse([],[],[],J+3*I,J)];
   C = [C1 ; C];                       
   D = [D1 ; D];    
   
   if   ~isnan(consIncluded)                        
  
        Mincluded = [consIncluded J+(1:K_Riv)];
        
        A = A(Mincluded,:);
        B = B(Mincluded,:);
        C = C(Mincluded,:);
        D = D(Mincluded,:);
           
   end
   
   
%%  Create matrix M for equality constraints: Mx = b

    if ~isnan(consIncluded)            
      
        M1 = [A C sparse([],[],[],length(consIncluded)+K_Riv,K_Flood)];    
        M2 = [-B sparse([],[],[],length(consIncluded)+K_Riv,2*J) ...
            sparse([],[],[],length(consIncluded)+K_Riv,K_Flood)];
        
    else
        
        M1 = [A C sparse([],[],[],J+K_Riv,K_Flood)];    
        M2 = [-B sparse([],[],[],J+K_Riv,2*J) ...
            sparse([],[],[],J+K_Riv,K_Flood)];
        
     end
    
        M = kron(speye(T_opt),M2) + kron(diag(ones(T_opt-1,1),-1),M1);    

%%  Create matrix M for inequality constraints: Nx <= d   
% Conduit flow < Storage, dt2*Qc - Sk <= 0 
% Irrigation flow < conduit flow, dt2*Qi - dt2*Qc <= 0 
% Water level - flood variable < Flood level
 
%     % System with irrigation inequality constraints, stage 1:
%     N1 = sparse([1:J 1:J J+(1:J) J+(1:J)],...
%         [1:J (2*J+K_Riv)+(1:J) (2*J+K_Riv)+(1:J) (3*J+K_Riv)+(1:J)],...
%         [-1*ones(1,J) ones(1,J) -ones(1,J) ones(1,J)],2*J,K);
%     N2 =  sparse([1:K_Flood  1:K_Flood ],...
%         [(J+3*seg_Flood(1:K_Flood)) (K+1 - (K_Flood:-1:1))],...
%         [ones(1,K_Flood) -ones(1,K_Flood)], K_Flood,K);
%     N = kron(speye(T_opt),[N1 ; N2]);
%     
    % System without irrigation inequality constraints, stage 2:
    N1 = sparse([1:J 1:J ],[1:J (2*J+K_Riv)+(1:J)],...
        [-1*ones(1,J) ones(1,J) ],J,K);
    N2 =  sparse([1:K_Flood  1:K_Flood ],...
        [(J+3*seg_Flood(1:K_Flood)) (K+1 - (K_Flood:-1:1))],...
        [ones(1,K_Flood) -ones(1,K_Flood)], K_Flood,K);

    N = kron(speye(T_opt),[N1 ; N2]);

    if  ~isnan(consIncluded)     

%         % System with irrigation inequality constraints, stage 1:
%         Nincluded = ((1:T_opt)-1)'*(2*J+K_Flood) + ...
%             [consIncluded J+consIncluded 2*J+(1:K_Flood) ];
        % System without irrigation inequality constraints, stage 2:
        Nincluded = ((1:T_opt)-1)'*(J+K_Flood) + ...
            [consIncluded J+(1:K_Flood) ];
        Nincluded = Nincluded';
        Nincluded = Nincluded(:);

        N = N(Nincluded,:);      
   
    end
       
    system.N = N;
    system.M = M;
    system.A = A;
    system.B = B;
    system.C = C;
    system.D = D;
    system.Q2 = Q2;
    system.J = J;
    system.I = I;
    system.K = K;
    system.K_Riv = K_Riv;
    
end
       


