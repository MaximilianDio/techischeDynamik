function [y, Dy, c, Dc] = subtask_e1(cs)
% simulation of crankshaft via manual partitioning of Constraint matrix and 
% analytical solution of beta

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
    %% define ODE for integration
    % define J and gamma
    J = [eye(1); -Ca\Cu];
    gamma = [0; -Ca\cs.dyn.ctt];
    
    M_hat = J'*cs.dyn.Mb*J;
    k_hat = J'*(cs.dyn.Mb*gamma + cs.dyn.kb);
    q_hat = J'*cs.dyn.qb;
    
    % ODE projected in valid movement space
    DDyu = matlabFunction(M_hat\(q_hat-k_hat),'vars',{[yu;ya;Dyu;Dya]});
    
    % x = [yu;Dyu]
    f = @(x,ya) [x(2); DDyu([x(1);ya(1);x(2);ya(2)])];
    
    %% solver options
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    
    %% iteration process
    for ii = 1:N
        %% apply boundary conditions -> solve for dependent coordinates
        y(ii,2) = asin(sin(y(ii,1))*cs.params.l1/cs.params.l2);
        if rank(Ca_(y(ii,:)')) == qc
            Dy(ii,2) = -Ca_(y(ii,:)')\(Cu_(y(ii,:)'))*Dy(ii,1);
        else
            error("dependent partition of constraint matrix is singular");
        end
        %% solve ODE
        if (ii+1 <= N)
%             [t,x] = BDFk(1,f,cs.sym.tspan(ii:ii+1),[y(ii,:), Dy(ii,:)]);
            [t,x] = ode45(@(t,x) f(x,[y(ii,2);Dy(ii,2)]),cs.sym.tspan(ii:ii+1),[y(ii,1); Dy(ii,1)],options);

            y(ii+1,1) = x(end,1);
            Dy(ii+1,1) = x(end,2);
        end
        
        c(ii,1) = c_f(y(ii,:)');
        Dc(ii,1) = Dc_f([y(ii,:),Dy(ii,:)]');
    end

end