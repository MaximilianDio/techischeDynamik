classdef CrankshaftTreeQR < CrankshaftTree
       
    properties
        name = "QR decomposition";
    end
    
    methods
        function obj = CrankshaftTreeQR(alpha0,Dalpha0)
            obj@CrankshaftTree(alpha0,Dalpha0);
        end    
    end
    methods 
        function results = solve(obj,tEnd)
            tic
            tspan = [0,tEnd];
            
            %%
            C = matlabFunction(obj.bc.C,'vars',{[obj.y;obj.Dy]});
            ct = matlabFunction(obj.bc.ct,'vars',{[obj.y;obj.Dy]});
            ctt = matlabFunction(obj.bc.ctt,'vars',{[obj.y;obj.Dy]});
            
            M = matlabFunction(obj.M,'vars',{[obj.y;obj.Dy]});
            k = matlabFunction(obj.k,'vars',{[obj.y;obj.Dy]});
            q = matlabFunction(obj.q,'vars',{[obj.y;obj.Dy]});
            
            %% solver options
            options = odeset('RelTol',obj.relTol,'AbsTol',obj.absTol);
            x0 = [obj.y0;obj.Dy0];
            
            %% solve as ODE
            [t,x] = ode45(@(t,x) obj.func(t,x,C,ct,ctt,M,k,q),tspan,x0,options);
            
            %% extract data from state
            
            % recalculate the dependent coordinates (pain in the ass but better than solving ODE in loop!)
            y =     x(:,1:2);
            Dy =    x(:,3:4);
            
            %% 
            disp(obj.name + " claculation time: " + string(toc));
            
            results = obj.validate(t,y,Dy);
        end
        function f = func(obj,t,x,C,ct,ctt,M,k,q)
            %% partitioning of dependent and independent coordiantes via QR
            % decomposition
            % C' = Q*R
            [Q,R] = qr(C(x)');

            % R = [R1;R2] (R2 = 0)
            R1 = R(1:size(R,2),:);
            % C' = Q1*R1 
            Q1 = Q(:,1:size(R,2));
            % J_ = Q2 ( projection matrix from yb to valid free movement of bunded system)
            Q2 = Q(:,size(R,2)+1:end);
            
            %% build other necessary matrices
            beta = @(x) 0;%-Q1*((R1')\ct(x)); % (no rehonom bonds)
            gamma = @(x) -Q1*((R1')\ctt(x));
            
            %% create ODE
            Mhat = @(x) Q2'*M(x)*Q2;
            khat = @(x) Q2'*(M(x)*gamma(x)+k(x));
            qhat = @(x) Q2'*q(x);

            Dyb = @(x) Q2*Q2'*x(3:4) + beta(x);
            DDyb = @(x) Q2*((Mhat(x))\(qhat(x)-khat(x)))+gamma(x);

            f = [Dyb(x) ;DDyb(x)];
            
        end
    end
end


