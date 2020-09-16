classdef CrankshaftTreeMPAnalytical < CrankshaftTree
       
    properties
        name = "MP - Analytical";
    end
    
    methods
        function obj = CrankshaftTreeMPAnalytical(alpha0,Dalpha0,solver)
            obj@CrankshaftTree(alpha0,Dalpha0,solver);
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
            DDyid = matlabFunction(Mhat\(qhat-khat),'vars',{[yid;yd;Dyid;Dyd]});
                        
            %% solver options
            options = odeset('RelTol',obj.relTol,'AbsTol',obj.absTol);
            x0 = [obj.y0(1);obj.Dy0(1)];
            
            %% solve as ODE
            [t,x] = obj.solver(@(t,x) obj.func(t,x,DDyid,Cid_,Cd_,qc),tspan,x0,options);
            
            Dx = zeros(size(x));
            for ii = 1:length(t)
                Dx(ii,:) = obj.func(t(ii),x(ii,:)',DDyid,Cid_,Cd_,qc)';
            end
            %% extract data from state
            ctt = matlabFunction(obj.bc.ctt,'vars',{[obj.y;obj.Dy]});
            % recalculate the dependent coordinates (pain in the ass but better than solving ODE in loop!)
            y =     [x(:,1), asin(sin(x(:,1))*obj.l1/obj.l2)];
            Dy =    [x(:,2), -(Cid_(y'))'./(Cd_(y'))'.*x(:,2)];
            DDy =   [Dx(:,2), -(Cid_(y'))'./(Cd_(y'))'.*Dx(:,2) - ctt([y';Dy'])'./(Cd_(y'))'];
            %%
            disp(obj.name + " claculation time: " + string(toc));
            
            results = obj.validate(t,y,Dy,DDy);
        end
        function f = func(obj,t,x,DDyi,Cid,Cd,qc)
            % x = [y;Dy];
            
            %% apply boundary conditions -> solve for dependent coordinates
                y = [x(1);
                     asin(sin(x(1))*obj.l1/obj.l2)];
                if rank(Cd(y)) == qc
                    Dy = [x(2); 
                          -Cd(y)\(Cid(y))*x(2)];
                else
                    error("dependent partition of constraint matrix is singular");
                end
            %% state space ode with correct dependent variables
            f = [Dy(1);
                DDyi([y;Dy])];
        end
    end
end


