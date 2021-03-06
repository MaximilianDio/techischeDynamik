classdef CrankshaftRedundantQR < CrankshaftRedundant
    
    properties 
        name = "QR - RC";
    end
    
    methods
        function obj = CrankshaftRedundantQR(alpha0,Dalpha0,solver)
            obj@CrankshaftRedundant(alpha0,Dalpha0,solver);
            
        end
        
        function results = solve(obj,tEnd)
            tic 
            
            tspan = 0:obj.tStep:tEnd;
            
            e = obj.e;
            q = obj.q;
            
            f = e-q;
            
            %%
            C = matlabFunction(obj.bc.C,'vars',{[obj.xI;obj.xII]});
            ct = @(x) obj.bc.ct;
            ctt = matlabFunction(obj.bc.ctt, 'vars',{[obj.xI;obj.xII]});
            
            M = @(x) obj.M;
            qc = @(x) obj.qc;%matlabFunction(qc, 'vars',{[xI;xII]});
            qe = @(x) obj.qe;%matlabFunction(qe, 'vars',{[xI;xII]});
            
            
            %% solver options
            options = odeset('RelTol',obj.relTol,'AbsTol',obj.absTol);
            x0 = [obj.xI0;obj.xII0];
            
            %% solve as ODE
            [t,x] = obj.solver(@(t,x) obj.func(t,x,C,ct,ctt,M,qc,qe),tspan,x0,options);
            
            %%
            disp(obj.name + " claculation time: " + string(toc));
            %% extract data from state
            Dx = zeros(size(x));
            for ii = 1:length(t)
                Dx(ii,:) = obj.func(t,x(ii,:)',C,ct,ctt,M,qc,qe)';
            end
            
            y = x(:,1:e);
            Dy = x(:,e+1:2*e);           
            DDy = Dx(:,e+1:2*e);
            
            results = obj.validate(t,y,Dy,DDy);
        end
        function Dy = func(obj,t,y,C,ct,ctt,M,qc,qe)
            % QR decomposition
            e = obj.e;
            q = obj.q;
            
            f = e-q;
            
            [Q,R] = qr(C(y)');

            Q1 = Q(:,1:end-f);
            Q2 = Q(:,end-f+1:end);

            R1 = R(1:end-f,:);

            beta = @(x) -Q1*((R1')\ct(x));
            gamma = @(x) -Q1*((R1')\ctt(x));

            Mhat = @(x) Q2'*M(x)*Q2;
            khat = @(x) Q2'*(M(x)*gamma(x)+qc(x));
            qhat = @(x) Q2'*qe(x);


            Dy = [Q2*Q2'*y(e+1:end)+beta(y);
                         Q2*(Mhat(y)\(qhat(y)-khat(y))) + gamma(y)];
        end
    end
    
end

