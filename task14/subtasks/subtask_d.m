function [y, Dy, c, Dc] = subtask_d(cs)
% simulation as DAE system
    syms lambda real;
    lambda0 = 0;
    
    DDy = [cs.dyn.Dyb;
                 cs.dyn.Mb\(cs.dyn.qb-cs.dyn.kb+cs.dyn.C'*lambda);
                 cs.dyn.c];

    
    f = matlabFunction(DDy,'vars',{[cs.dyn.yb;cs.dyn.Dyb;lambda]});
    
    % DAE by choosing last entry to have no dynamics
    M = diag([1,1,1,1,0]);
    options = odeset('Mass',M,'RelTol',1e-4,'AbsTol',1e-6);
    

    y0 = [cs.sym.alpha0;cs.sym.beta0;cs.sym.Dalpha0;cs.sym.Dbeta0;lambda0];
    [t,y] = ode23t(@(t,y) f(y),cs.sym.tspan,y0,options);
    
    y = 0;
    Dy = 0;
    c = 0;
    Dc = 0;
    % differential algebraic equation:
    % Dy_b = z_b
    % M_b*Dz_b + k_b - C'*lambda =  q_b
    % c = 0
%     Dx = @(t,x) []; 
    
%     [t,y] = ode15s(Dy,tspan,x0);

end