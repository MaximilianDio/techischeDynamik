classdef CrankshaftTree < CrankshaftClass
    
    properties
        % coordinates in tree structure
        y; 
        Dy;
        
        % initial condition
        y0;
        Dy0;
        
        % matrices for principle of D'Alembert
        M;
        k;
        q;
        
    end
    
    methods
        function obj = CrankshaftTree(alpha0,Dalpha0)
            obj@CrankshaftClass(alpha0,Dalpha0);
            
            y = sym('y',[2,1],'real');
            Dy = sym('Dy',[2,1],'real');
            
            obj.y = y;
            obj.Dy = Dy;
            
            %% define initial conditions
            obj.y0 = [obj.alpha0;obj.beta0];
            obj.Dy0 = [obj.Dalpha0;obj.Dbeta0];
            
            %% define required matrices and vectors for principle of D'Alembert
            % Mass matrix
            obj.M = [(obj.m1+obj.m2)*obj.l1^2               -obj.m2*obj.l1*obj.l2*cos(y(1)+y(2));  
                     -obj.m2*obj.l1*obj.l2*cos(y(1)+y(2))   obj.m2*obj.l2^2];
            % vector of generalized coriolis, centrifugal and gyroscopic forces
            obj.k = [obj.m2*obj.l1*obj.l2*Dy(2)^2*sin(y(1)+y(2)); 
                     obj.m2*obj.l1*obj.l2*Dy(1)^2*sin(y(1)+y(2))]; 
            % vector of generalized embedded forces
            obj.q = [-(obj.m1+obj.m2)*obj.g*obj.l1*cos(y(1));
                     obj.m2*obj.g*obj.l2*cos(y(2))];
                       
        end
        function [xI_1, xI_2] = position(obj,y)
            xI_1 = [obj.l1*cos(y(:,1)),...
                   obj.l1*sin(y(:,1))];
            xI_2 = [obj.l1*cos(y(:,1))+obj.l2*cos(y(:,2)),...
                   obj.l1*sin(y(:,1))-obj.l2*sin(y(:,2))];
        end
        function [xII_1, xII_2] = velocity(obj,y,Dy)
            xII_1 = [-obj.l1*sin(y(:,1)).*Dy(:,1),...
                     obj.l1*cos(y(:,1)).*Dy(:,1)];
            xII_2 = [-obj.l1*sin(y(:,1)).*Dy(:,1) - obj.l2*sin(y(:,2)).*Dy(:,2),...
                     obj.l1*cos(y(:,1)).*Dy(:,1)-obj.l2*cos(y(:,2)).*Dy(:,2)];
        end
        function [alpha,beta] = angles(obj,y)
%             toPi = @(angle) angle - 2*pi*floor( (angle+pi)/(2*pi) );
            alpha = wrapToPi(y(:,1)+pi/2)-pi/2;
            beta = wrapToPi(y(:,2)+pi/2)-pi/2;
        end
        
        function [Dalpha,Dbeta] = angleVelocities(obj,Dy)
            Dalpha = Dy(:,1);
            Dbeta = Dy(:,2);
        end
        
    end
    methods (Access = protected)
        function results = validate(obj,t,y,Dy)
            results.info.name = obj.name;
            % 
            results.t = t;
            results.y = y;
            results.Dy = Dy;
            
            % validate boundary conditions
            c_f = matlabFunction(obj.bc.c,'vars',{[obj.y;obj.Dy]});
            Dc_f = matlabFunction(obj.bc.Dc,'vars',{[obj.y;obj.Dy]});
            
            T = length(t);
            
            c = zeros(T,1);
            Dc = zeros(T,1);
            for ii = 1:T        
                c(ii,1) = c_f([y(ii,:),Dy(ii,:)]');
                Dc(ii,1) = Dc_f([y(ii,:),Dy(ii,:)]');
            end
            
            results.c = c;
            results.Dc = Dc;
        end
    end
end


