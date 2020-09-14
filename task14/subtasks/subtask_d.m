function [x, Dy, c, Dc, task_info] = subtask_d(cs)
% simulation as DAE system

    %% Task information
    task_info.name = "DAE desctription";
    
   	%%
    N = length(cs.sym.tspan);
    
    % pre-allocate vectors
    y = zeros(N,2);
    Dy = zeros(N,2);
    c = zeros(N,1);
    Dc = zeros(N,1);
    
    % fill in independent initial conditions
    y(1,:) = [cs.sym.alpha0;cs.sym.beta0];
    Dy(1,:) = [cs.sym.Dalpha0;cs.sym.Dbeta0];
    %% function handle to evaluate error to contraint
    c_f = matlabFunction(cs.dyn.c,'vars',{[cs.dyn.yb]});
    Dc_f = matlabFunction(cs.dyn.Dc,'vars',{[cs.dyn.yb;cs.dyn.Dyb]});
    
    %%
    syms t lambda real;
    lambda0 = 0;
    
    % define DAE and singular mass matrix s.t. DAE has Index 1
    f = [cs.dyn.Dyb;
                 cs.dyn.Mb\(cs.dyn.qb-cs.dyn.kb+cs.dyn.C'*lambda);
                 -cs.dyn.ctt];

    M = [eye(4),zeros(4,1); 0,0 cs.dyn.C,0];
    
    f = matlabFunction(f,'vars',{t,[cs.dyn.yb;cs.dyn.Dyb;lambda]});
    M = matlabFunction(M,'vars',{t,[cs.dyn.yb;cs.dyn.Dyb;lambda]});
    
    
    options = odeset('Mass',M,'RelTol',1e-10,'AbsTol',1e-10);
    

    x0 = [y(1,:),Dy(1,:),lambda0];
    [~,x] = ode15s(f,cs.sym.tspan,x0,options);
    
    y = x(:,1:2);
    Dy = x(:,3:4);
    
    for ii = 1:N        
        c(ii,1) = c_f(y(ii,:)');
        Dc(ii,1) = Dc_f([y(ii,:),Dy(ii,:)]');
    end
    


end