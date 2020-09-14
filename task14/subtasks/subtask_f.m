function [y, Dy, c, Dc, task_info] = subtask_f(cs)
%%  subtask e2: simulation via manual coordinate partitioning by integration
    tspan = cs.sym.tspan;
%     tspan = cs.sym.tspan(1):0.0005:cs.sym.tspan(end);
    %% Task information
    task_info.name = "QR-Decomposition";
    
    %% preprcessing
    N = length(tspan);
    
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
    
    %% convert symbolic functions to function handles    
    C = matlabFunction(cs.dyn.C,'vars',{[cs.dyn.yb]});
    ctt = matlabFunction(cs.dyn.ctt,'vars',{[cs.dyn.yb;cs.dyn.Dyb]});
    
    Mb = matlabFunction(cs.dyn.Mb,'vars',{[cs.dyn.yb;cs.dyn.Dyb]});
    kb = matlabFunction(cs.dyn.kb,'vars',{[cs.dyn.yb;cs.dyn.Dyb]});
    qb = matlabFunction(cs.dyn.qb,'vars',{[cs.dyn.yb;cs.dyn.Dyb]});
    
    %% solver options
    options = odeset('RelTol',1e-8,'AbsTol',1e-8);
    
    for ii = 1:N
        %% partitioning of dependent and independent coordiantes via QR
    	% decomposition
        % C' = Q*R
        [Q,R] = qr(C(y(ii,:)')');
        
        % R = [R1;R2] (R2 = 0)
        R1 = R(1:size(R,2),:);
        % C' = Q1*R1 
        Q1 = Q(:,1:size(R,2));
        % J_ = Q2 ( projection matrix from yb to valid free movement of bunded system)
        Q2 = Q(:,size(R,2)+1:end);
        
        %% build other necessary matrices
        beta = @(x) 0 ;% (no rehonom bonds)
        gamma = @(x) -Q1*((R1')\ctt(x));
        
        %% create ODE
        M_hat = @(x) Q2'*Mb(x)*Q2;
        k_hat = @(x) Q2'*(Mb(x)*gamma(x)+kb(x));
        q_hat = @(x) Q2'*qb(x);
        
        Dyb = @(x) Q2*Q2'*[x(3);x(4)] + beta(x);
        DDyb = @(x) Q2*((M_hat(x))\(q_hat(x)-k_hat(x)))+gamma(x);
        
        DDx = @(x) [Dyb(x) ;DDyb(x)];
        
        %% solve ODE
        if ii+1 <= N
            [~,x] = ode15s(@(t,x) DDx(x),tspan(ii:ii+1),[y(ii,:), Dy(ii,:)]',options);

            y(ii+1,:) = x(end,1:2);
            Dy(ii+1,:) = x(end,3:4);

        end
        
        c(ii,1) = c_f(y(ii,:)');
        Dc(ii,1) = Dc_f([y(ii,:),Dy(ii,:)]');
    end
    
    y = interp1(tspan,y,cs.sym.tspan);
    Dy = interp1(tspan,Dy,cs.sym.tspan);
    c = interp1(tspan,c,cs.sym.tspan);
    Dc = interp1(tspan,Dc,cs.sym.tspan);
    
    

end