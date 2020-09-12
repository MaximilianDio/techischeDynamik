function [y, Dy, c, Dc, task_info] = subtask_i(cs)
%%  subtask e2: simulation via manual coordinate partitioning by integration
    
    %% Task information
    task_info.name = "subtask i";
    
    %% preprcessing
    N = length(cs.sym.tspan);
    
    % pre-allocate vectors
    y = zeros(N,4);
    Dy = zeros(N,4);
    c = zeros(N,1);
    Dc = zeros(N,1);

    %% build DAE
    M = cs.DAE.M;
    qe = cs.DAE.qe;
    C =  matlabFunction(cs.DAE.C,'vars',{cs.DAE.x});
    
    q = 3; % number of constraints
    p = 2; % number of bodies
    
    f = @(t,x) [x(p*2+1:p*2*2);M\(qe+C(x(1:p*2))'*x(p*2*2+1:p*2*2+q));zeros(q,1)];
    
    MM = @(t,x) [eye(p*2*2,p*2*2+q);zeros(q,p*2),C(x(1:p*2)),zeros(q)];
    
   
    
    %% solve DAE
    % initial values:
    lambda0 = [0 0 0]; %reaction forces
      
    % solve DAE
    options = odeset('Mass',MM,'RelTol',1e0,'AbsTol',1e0);
    

    x0 = [cs.DAE.x0 ,cs.DAE.Dx0,lambda0]';
    [~,x] = ode15s(f,cs.sym.tspan,x0,options);
    
    %% postprecessing
    
end