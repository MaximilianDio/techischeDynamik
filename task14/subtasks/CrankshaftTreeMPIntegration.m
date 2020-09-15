classdef CrankshaftTreeMPIntegration < CrankshaftTree
       
    properties
        name = "MP - Integration";
    end
    
    methods
        function obj = CrankshaftTreeMPIntegration(alpha0,Dalpha0)
            obj@CrankshaftTree(alpha0,Dalpha0);
        end    
    end
    methods 
        function results = solve(obj,tEnd)
            tic
            
            tspan = 0:obj.tStep:tEnd;
            
            % independent variable (alpha)
            yid = obj.y(1);
            Dyid = obj.Dy(1);
            % dependent variable (beta)
            yd = obj.y(2);
            Dyd = obj.Dy(2);
            
            % partitioning of C in independent and dependent matrices
            Cid = obj.bc.C(1);
            Cd = obj.bc.C(2);
            Cid_ = matlabFunction(Cid,'vars',{[obj.y]});
            Cd_ = matlabFunction(Cd,'vars',{[obj.y]});
            
            % number of constraints
            qc = length(yd);
            %% define projections and projected matrices
            % define J and gamma
            J = [eye(1); -Cd\Cid];
            gamma = [0; -Cd\obj.bc.ctt];

            Mhat = J'*obj.M*J;
            khat = J'*(obj.M*gamma + obj.k);
            qhat = J'*obj.q;
            
            %% ODE projected in independent coordinate space
            DDy = matlabFunction([Dyid;
                                  -(Cd\Cid)*Dyid;
                                  Mhat\(qhat-khat);
                                  0],'vars',{[yid;yd;Dyid;Dyd]});
                        
            %% solver options
            options = odeset('RelTol',obj.relTol,'AbsTol',obj.absTol);
            x0 = [obj.y0;obj.Dy0];
            
            %% solve as ODE
            [t,x] = ode45(@(t,x) obj.func(t,x,DDy,Cid_,Cd_,qc),tspan,x0,options);
            
            %% extract data from state
            
            % recalculate the dependent coordinates (pain in the ass but better than solving ODE in loop!)
            y =    x(:,1:2);
            Dy =   [x(:,3), -(Cid_(y'))'./(Cd_(y'))'.*x(:,3)];
            
            %%
            disp(obj.name + " claculation time: " + string(toc));
            
            results = obj.validate(t,y,Dy);
        end
        function f = func(obj,t,x,DDy,Cid,Cd,qc)
            % x = [y;Dy];
            
            %% apply boundary conditions -> solve for dependent coordinates
                if rank(Cd(x(1:2))) == qc
                   x(4) = -Cd(x(1:2))\(Cid(x(1:2)))*x(3);
                else
                    error("dependent partition of constraint matrix is singular");
                end
            %% state space ode with correct dependent variables
            f = DDy(x);
        end
        
    end
end


