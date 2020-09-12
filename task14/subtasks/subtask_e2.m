function [y, Dy, c, Dc, task_info] = subtask_e2(cs)
%%  subtask e2: simulation via manual coordinate partitioning by integration
    
    %% Task information
    task_info.name = "subtask e2";
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
    
    %%
    % function to evaluate error to contraint
    c_f = matlabFunction(cs.dyn.c,'vars',{[cs.dyn.yb]});
    Dc_f = matlabFunction(cs.dyn.Dc,'vars',{[cs.dyn.yb;cs.dyn.Dyb]});
    
    %%
    % partitioning of dependent and independent coordiantes
    % alpha is taken as independent coordinate:
    
    % independent variable (alpha)
    yu = cs.dyn.yb(1);
    Dyu = cs.dyn.Dyb(1);
    % dependent variable (beta)
    ya = cs.dyn.yb(2);
    Dya = cs.dyn.Dyb(2);
    
    Cu = cs.dyn.C(1);
    Ca = cs.dyn.C(2);
    Cu_ = matlabFunction(Cu,'vars',{[cs.dyn.yb]});
    Ca_ = matlabFunction(Ca,'vars',{[cs.dyn.yb]});
    
    qc = length(ya);
    %% define projections and projected matrices
    % define J and gamma
    J = [eye(1); -Ca\Cu];
    gamma = [0; -Ca\cs.dyn.ctt];
    
    M_hat = J'*cs.dyn.Mb*J;
    k_hat = J'*(cs.dyn.Mb*gamma + cs.dyn.kb);
    q_hat = J'*cs.dyn.qb;
    
    %% ODE projected in valid movement space
    
    DDyb = matlabFunction([Dyu; 
                           -(Ca\Cu)*Dyu;     
                           M_hat\(q_hat-k_hat);0],'vars',{[yu;ya;Dyu;Dya]});
    
    %% iteration process - numerical solution of ODE
    for ii = 1:N
        if rank(Ca_(y(ii,:)')) == qc
            Dy(ii,2) = -Ca_(y(ii,:)')\(Cu_(y(ii,:)'))*Dy(ii,1);
        else
            error("dependent partition of constraint matrix is singular");
        end
        
        if ii+1 <= N
            [~,x] = ode45(@(t,x) DDyb(x),cs.sym.tspan(ii:ii+1),[y(ii,:), Dy(ii,:)]');

            y(ii+1,:) = x(end,1:2);
            Dy(ii+1,1) = x(end,3);

        end
        
        c(ii,1) = c_f(y(ii,:)');
        Dc(ii,1) = Dc_f([y(ii,:),Dy(ii,:)]');
    end
    

end