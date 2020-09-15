classdef Crankshaft<handle
    %CRANKSHAFTCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties        
        % initial conditions
        alpha0;         %rad
        beta0; 
        Dalpha0;        %rad/s
        Dbeta0;
        
        %boundary conditions
        bc;
    end
    properties (Constant)
       % point masses
        m1 = .1;        %kg
        m2 =  1.0;      %kg
        % rod lengths
        l1 = 0.3;       %m
        l2 = 1;         %m
        % gravitational constant
        g = 9.81;       % m/s^2
    
        % simulation
        relTol = 1e-4;
        absTol = 1e-4;
        tStep = 0.01; 
    end
    properties (Abstract)
        name;
    end
    
    methods 
        function obj = Crankshaft(alpha0,Dalpha0)
            obj.alpha0 = alpha0;
            obj.Dalpha0 = Dalpha0;
            
            %% calculate Initial conditions in form aof agles
            y = sym('y',[2,1],'real');
            Dy = sym('Dy',[2,1],'real');
            
            %% generate boundary condition
            c = sin(y(1))*obj.l1-sin(y(2))*obj.l2;
            obj.bc = BoundaryCondition(c,y,Dy);
            
            %% solve initial values
            tmp = solve(subs(obj.bc.c,y(1),obj.alpha0) == 0,y(2));
            obj.beta0 = double(tmp(1));
            obj.Dbeta0 = double(solve(subs(obj.bc.Dc,{y(1),y(2),Dy(1)},{alpha0,obj.beta0,Dalpha0}) == 0,Dy(2)));
            
        end
    end
    methods (Abstract)
        results = solve(obj,tEnd);
        [xI_1, xI_2] = position(obj,y);
        [xII_1, xII_2] = velocity(obj,y,Dy);
        [DxII_1, DxII_2] = acceleration(obj,y,Dy);
        [alpha,beta] = angles(obj,y);
        [Dalpha,Dbeta] = angleVelocities(obj,Dy);
    end
end

