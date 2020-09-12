function [y, Dy, c, Dc, task_info] = subtask_f(cs)
%%  subtask e2: simulation via manual coordinate partitioning by integration
    
    %% Task information
    task_info.name = "subtask f";
    
    %% preprcessing
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
    
    %% convert symbolic functions to function handles    
    C = matlabFunction(cs.dyn.C,'vars',{[cs.dyn.yb]});
    ctt = matlabFunction(cs.dyn.ctt,'vars',{[cs.dyn.yb;cs.dyn.Dyb]});
    
    Mb = matlabFunction(cs.dyn.Mb,'vars',{[cs.dyn.yb;cs.dyn.Dyb]});
    kb = matlabFunction(cs.dyn.kb,'vars',{[cs.dyn.yb;cs.dyn.Dyb]});
    qb = matlabFunction(cs.dyn.qb,'vars',{[cs.dyn.yb;cs.dyn.Dyb]});
    

    
    for ii = 1:N
        %% partitioning of dependent and independent coordiantes via QR
    	% decomposition
        % C' = Q*R
        [Q,R] = qr(C(y(ii,:)')');
        
        % R = [R1;R2] (R2 = 0)
        R1 = R(1:size(R,2),:);
        % C' = Q1*R1 
        Q1 = Q(:,1:size(R,2));
        % J_ = Q2 ( projection matrix from yb to valid movement of bunded system)
        J_ = Q(:,size(R,2)+1:end);
        
        %% build other necessary matrices
        beta = @(x) 0 ;% (no rehonom bonds)
        gamma = @(x) -Q1*(R1'\ctt(x));
        
        %% create ODE
        M_hat = @(x) J_'*Mb(x)*J_;
        k_hat = @(x) J_'*(Mb(x)*gamma(x)+kb(x));
        q_hat = @(x) J_'*qb(x);
        
        Dyb = @(x) J_*J_'*[x(1);x(2)] + beta(x);
        DDyb = @(x) J_*(M_hat(x)\(q_hat(x)-k_hat(x)))+gamma(x);
        
        DDx = @(x) [Dyb(x) ;DDyb(x)];
        
        %% solve ODE
        if ii+1 <= N
%             h = cs.sym.tspan(ii+1)-cs.sym.tspan(ii);
%             x = ([y(ii,:), Dy(ii,:)]' + h*1/2*DDx([y(ii,:), Dy(ii,:)]'+h/2*DDx([y(ii,:), Dy(ii,:)]')))';
            [~,x] = ode45(@(t,x) DDx(x),cs.sym.tspan(ii:ii+1),[y(ii,:), Dy(ii,:)]');

            y(ii+1,:) = x(end,1:2);
            Dy(ii+1,:) = x(end,3:4);

        end
        
        c(ii,1) = c_f(y(ii,:)');
        Dc(ii,1) = Dc_f([y(ii,:),Dy(ii,:)]');
    end
    
    

end