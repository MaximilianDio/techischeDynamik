classdef CrankshaftRedundantDAE < CrankshaftRedundant
    
    properties 
        name = "DAE - RC";
    end
    
    methods
        function obj = CrankshaftRedundantDAE(alpha0,Dalpha0,solver)
            obj@CrankshaftRedundant(alpha0,Dalpha0,solver);
            
        end
        
        function results = solve(obj,tEnd)
            tic 
            
            tspan = 0:obj.tStep:tEnd;
            
            e = obj.e;
            q = obj.q;

            %%
            syms  F1 F2 F3 F4 F5 real
            lambda = [F1; F2; F3; F4; F5];

            %% DAE formulation
            DxI = matlabFunction(obj.Z*obj.xII, 'vars',{[obj.xI;obj.xII;lambda]});
            DxII = matlabFunction(obj.M\(obj.qe-obj.qc+obj.bc.C'*lambda),'vars',{[obj.xI;obj.xII;lambda]});
            ctt = matlabFunction(obj.bc.ctt, 'vars',{[obj.xI;obj.xII;lambda]});
            C = matlabFunction(obj.bc.C,'vars',{[obj.xI;obj.xII;lambda]});
            M_ = @(t,x) [eye(e) zeros(e) zeros(e,q); 
                       zeros(e) eye(e) zeros(e,q); 
                       zeros(q,e) C(x)    zeros(q)];

            f = @(t,x) [DxI(x);DxII(x);-ctt(x)];
            
            
            %%
            lambda0 = zeros(size(lambda));
            x0 = [obj.xI0;obj.xII0;lambda0];
            
            options = odeset("Mass",M_,"RelTol",obj.relTol,"AbsTol",obj.absTol);            
            
            [t,x] = obj.solver(f,tspan,x0,options);
            %%
            disp(obj.name + " claculation time: " + string(toc));
            
            %% extract data from state
            Dx = zeros(size(x));
            for ii = 1:length(t)
                Dx(ii,:) = f(t(ii),x(ii,:)')';
            end
            
            y = x(:,1:e);
            Dy = x(:,e+1:2*e);
            Dyy = Dx(:,e+1:2*e);
       
            % process 
            results = obj.validate(t,y,Dy,Dyy);
        end
    end
    
end

