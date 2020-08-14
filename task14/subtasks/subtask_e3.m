function [y, Dy, c, Dc] = subtask_e2(cs)
    N = length(cs.sym.tspan);
    
    % pre-allocate vectors
    y = zeros(N,2);
    Dy = zeros(N,2);
    c = zeros(N,1);
    Dc = zeros(N,1);
    % fill in independent initial conditions
    y(1,1) = cs.sym.alpha0;
    Dy(1,1) = cs.sym.Dalpha0;
    % step size
    h = cs.sym.tspan(2)-cs.sym.tspan(1);
    
    %%
    % function to evaluate error to contraint
    c_f = matlabFunction(cs.dyn.c,'vars',{[cs.dyn.yb]});
    Dc_f = matlabFunction(cs.dyn.Dc,'vars',{[cs.dyn.yb;cs.dyn.Dyb]});
    
    %%
    % partitioning of dependent and independent coordiantes
    % alpha is taken as independent coordinate:
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
    
    
    DDy = matlabFunction(M_hat\(q_hat-k_hat),'vars',{[cs.dyn.yb(1);cs.dyn.Dyb(1)],[cs.dyn.yb(2);cs.dyn.Dyb(2)]});
    f = @(t,y,ya) [y(2);DDy(y,ya)];
    %     f = @(t,x) MKSode(t,x,DDy);
    
    
    options = odeset('RelTol',1e-5,'AbsTol',1e-5);
    %% iteration process
    for ii = 1:N
        %% apply boundary conditions
        y(ii,2) = asin(sin(y(ii,1))*cs.params.l1/cs.params.l2);
        try
            Dy(ii,2) = -Ca_(y(ii,:)')\(Cu_(y(ii,:)'))*Dy(ii,1);
        catch
            error("dependent partition of constraint matrix is singular");
        end
        %% solve ODE
        f_ = @(t,x) f(t,x,[y(ii,2);Dy(ii,2)]);
        if (ii+1 <= N)
%             [t,x] = BDFk(1,f,cs.sym.tspan(ii:ii+1),[y(ii,:), Dy(ii,:)]);
            [t,x] = ode45(f_,cs.sym.tspan(ii:ii+1),[y(ii,1), Dy(ii,1)],options);

            y(ii+1,:) = x(end,:);
        end
        
        c(ii,1) = c_f(y(ii,:)');
        Dc(ii,1) = Dc_f([y(ii,:),Dy(ii,:)]');
    end
end