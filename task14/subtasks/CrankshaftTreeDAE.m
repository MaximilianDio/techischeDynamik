classdef CrankshaftTreeDAE < CrankshaftTree
       
    properties
        name = "DAE - tree coordinates";
    end
    
    methods
        function obj = CrankshaftTreeDAE(alpha0,Dalpha0)
            obj@CrankshaftTree(alpha0,Dalpha0);
        end    
    end
    methods 
        function results = solve(obj,tEnd)
            tic
            
            tspan = 0:obj.tStep:tEnd;
            
            syms t lambda real;
            lambda0 = 0;
            
            %% DAE with reduced index 
            % define right handside of DAE
            f = [obj.Dy;
                        obj.M\(obj.q-obj.k+obj.bc.C'*lambda);
                        -obj.bc.ctt];
            % define left handside of DAE (mass matrix)
            M = [eye(4),zeros(4,1); 0,0 obj.bc.C,0];
            
            % convert to function handle
            f = matlabFunction(f,'vars',{t,[obj.y;obj.Dy;lambda]});
            M = matlabFunction(M,'vars',{t,[obj.y;obj.Dy;lambda]});
            
            %% solver options
            options = odeset('Mass',M,'RelTol',obj.relTol,'AbsTol',obj.absTol);
            x0 = [obj.y0;obj.Dy0;lambda0];
            
            %% solve as DAE
            [t,x] = ode15s(f,tspan,x0,options);
            
            Dx = f(t,x')';
            
            %% extract data from state
            y = x(:,1:2);
            Dy = x(:,3:4);
            DDy = Dx(:,3:4);
            
            %%
            disp(obj.name + " claculation time: " + string(toc));
            
            results = obj.validate(t,y,Dy,DDy);
        end
        
    end
end

