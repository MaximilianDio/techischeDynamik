%% define matrices for Analysis in Tree structure of crankshaft
function crankshaft = crankshaft(alpha,beta,Dalpha,Dbeta)
    % variables for convenience
	% point masses
    m1 = .1;       %kg
    m2 =  1;       %kg
    % rod lengths
    l1 = 0.3;       %m
    l2 = 1;         %m
    % gravitational constant
    g = 9.81;           % m/s^2
    
    % initial conditions
    alpha0 = 0.1;   %rad
    Dalpha0 = 0.1;  %rad/s
    tspan = 0:0.01:2;
      
    %% system parameters 
        crankshaft.params.m1 = m1;
        crankshaft.params.m2 = m2;
        crankshaft.params.l1 = l1;
        crankshaft.params.l2 = l2;
    
    %% environment parameters
        crankshaft.env.g = g;
    
    %% kinematics
        crankshaft.kin.r1 = [l1*cos(alpha),l1*sin(alpha)];
        crankshaft.kin.r2 = [l1*cos(alpha)+l2*cos(beta),l1*sin(alpha)-l2*sin(beta)];
        
        crankshaft.kin.v1 = [-l1*sin(alpha)*Dalpha,l1*cos(alpha)*Dalpha];
        crankshaft.kin.v2 = [-l1*sin(alpha)*Dalpha - l2*sin(beta)*Dbeta,...
                                l1*cos(alpha)*Dalpha-l2*cos(beta)*Dbeta];
        
    %% dynamics        
        % general coordinates in tree structure 
        % NOTE: beta and alpha are not independent when crankshaft has to
        % enforce kinematic boundary condition
        crankshaft.dyn.yb = [alpha; beta];
        crankshaft.dyn.Dyb = [Dalpha; Dbeta];
        % mass matrix
        crankshaft.dyn.Mb = [(m1+m2)*l1^2                   -m2*l1*l2*cos(alpha+beta);  
                             -m2*l1*l2*cos(alpha+beta)       m2*l2^2];
        % vector of generalized coriolis, centrifugal and gyroscopic forces
        crankshaft.dyn.kb = [m2*l1*l2*Dbeta^2*sin(alpha+beta); 
                             m2*l1*l2*Dalpha^2*sin(alpha+beta)]; 
        % vector of generalized embedded forces
        crankshaft.dyn.qb = [-(m1+m2)*g*l1*cos(alpha);
                                  m2*g*l2*cos(beta)];
        
	%% boundary conditions
        % kinematic boundary condition (location)
        crankshaft.dyn.c = sin(alpha)*l1-sin(beta)*l2; % =0
        % jacobian matrix of c dc/dy
        crankshaft.dyn.C = jacobian(crankshaft.dyn.c);
        
        % kinematic boundary condition (velocity)
        crankshaft.dyn.Dc = crankshaft.dyn.C*[Dalpha;Dbeta]; % =0
        
        crankshaft.dyn.ct = 0;
        crankshaft.dyn.ctt = [-sin(alpha)*l1 sin(beta)*l2]*[Dalpha^2;Dbeta^2];
    
    %% simulation parameters
        % initial conditions 
        crankshaft.sym.alpha0 = alpha0;
        crankshaft.sym.Dalpha0 = Dalpha0;
        % solve for initial conditions of beta
        tmp = solve(subs(crankshaft.dyn.c,'alpha',alpha0) == 0,beta);
        crankshaft.sym.beta0 = double(tmp(1));
        crankshaft.sym.Dbeta0 = double(solve(subs(crankshaft.dyn.Dc,{'alpha','beta','Dalpha'},{alpha0,crankshaft.sym.beta0,Dalpha0}) == 0,Dbeta));

        % time span
        crankshaft.sym.tspan = tspan;
        
    %% DAE formulation (for subtask i and j)
        syms x1 y1 x2 y2 real
        
        crankshaft.DAE.x = [x1; y1; x2; y2];
        
        crankshaft.DAE.x0 = double([subs(crankshaft.kin.r1,{'alpha','beta'},...
            [crankshaft.sym.alpha0,crankshaft.sym.beta0]),...
            subs(crankshaft.kin.r2,{'alpha','beta'},...
            [crankshaft.sym.alpha0,crankshaft.sym.beta0])]);
        
        crankshaft.DAE.Dx0 = double([subs(crankshaft.kin.v1,{'alpha','beta','Dalpha','Dbeta'},...
            [crankshaft.sym.alpha0,crankshaft.sym.beta0,crankshaft.sym.Dalpha0,crankshaft.sym.Dbeta0]),...
            subs(crankshaft.kin.v2,{'alpha','beta','Dalpha','Dbeta'},...
            [crankshaft.sym.alpha0,crankshaft.sym.beta0,crankshaft.sym.Dalpha0,crankshaft.sym.Dbeta0])...
            ]);
        
        crankshaft.DAE.M = diag([m1,m1,m2,m2]);
        crankshaft.DAE.qe = [0;-m1*g;0;-m2*g];
        
        crankshaft.DAE.c = [x1^2+y1^2-l1;(x1-x2)^2+(y1-y2)^2-l2;y2];
        crankshaft.DAE.C = [2*x1        2*y1        0           0;
                            2*(x1-x2)   2*(y1-y2)   -2*(x1-x2)  -2*(y1-y2);
                            0           0           0           1];

end