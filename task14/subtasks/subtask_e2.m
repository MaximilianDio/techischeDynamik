function [y, Dy, c, Dc] = subtask_e2(cs)
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
    
    
    
    %% define ODE for integration
    % define J and gamma
    J = [eye(1); -Ca\Cu];
    gamma = [0; -Ca\cs.dyn.ctt];
    
    M_hat = J'*cs.dyn.Mb*J;
    k_hat = J'*(cs.dyn.Mb*gamma + cs.dyn.kb);
    q_hat = J'*cs.dyn.qb;
    
    % ODE in 
    DDyu = M_hat\(q_hat-k_hat);
    
    DDx = matlabFunction([Dyu;Dya;DDyu;-Ca\Cu*Dyu-Ca\cs.dyn.ctt],'vars',{[yu;ya;Dyu;Dya]});
    f = @(t,x) DDx(x);    
    
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [t,x] = ode45(f,cs.sym.tspan,[y(1,:), Dy(1,:)]',options);
    
    y = x(:,1:2);
    Dy = x(:,3:4);
    %% iteration process
%     for ii = 1:N
%         
%         c(ii,1) = c_f(y(ii,:)');
%         Dc(ii,1) = Dc_f([y(ii,:),Dy(ii,:)]');
%     end
end