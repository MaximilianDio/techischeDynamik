classdef CrankshaftRedundant < CrankshaftClass
    
    properties
        % coordinates in tree structure
        xI; 
        xII;
        
        % initial condition
        xI0;
        xII0;
        
        % matrices for principle of D'Alembert
        M;
        qe;
        qc;
        
        Z;
        
        e; % degrees of freedom unconstrained
        q; % number of total constaints
        
    end
    
    methods
        function obj = CrankshaftRedundant(alpha0,Dalpha0)
            obj@CrankshaftClass(alpha0,Dalpha0);
            
            syms x1 y1 z1 Dx1 Dy1 Dy1 Dz1 x2 y2 z2 Dx2 Dy2 Dy2 Dz2 real
 
            xI = [x1; y1; z1; x2; y2; z2];
            xII = [Dx1; Dy1; Dz1; Dx2; Dy2; Dz2];
            
            obj.xI = xI;
            obj.xII = xII;
            
            p = 2;
            ei = 3; 
            qi = 2;

            obj.e = ei*p;
            obj.q = qi*p+1;
            
            obj.Z = eye(obj.e);
            
            %% define required matrices and vectors for principle of D'Alembert
            obj.M = blkdiag(eye(ei)*obj.m1,eye(ei)*obj.m2);
            obj.qc = [0;0;0;0;0;0];
            obj.qe = [0; 0; -obj.m1*obj.g;0; 0; -obj.m2*obj.g];
            
            %% define initial conditions
            obj.xI0 =  [ 0; 
                         obj.l1*cos(alpha0);
                        obj.l1*sin(alpha0);
                        0;
                        obj.l1*cos(alpha0)+obj.l2*cos(obj.beta0);
                        obj.l1*sin(alpha0)-obj.l2*sin(obj.beta0)];
                obj.xII0 = [ 0;
                            -obj.l1*sin(alpha0).*Dalpha0;
                            obj.l1*cos(alpha0).*Dalpha0;
                            0;
                            -obj.l1*sin(alpha0).*Dalpha0 - obj.l2*sin(obj.beta0).*obj.Dbeta0;
                            obj.l1*cos(alpha0).*Dalpha0-obj.l2*cos(obj.beta0).*obj.Dbeta0];
            
            %% obj.
            obj.bc.c = [x1;
                z1^2+y1^2-obj.l1^2;
                x2;
                (y2-y1)^2+(z2-z1)^2-obj.l2^2;
                z2];
            obj.bc.C = jacobian(obj.bc.c,xI);
            obj.bc.Dc = obj.bc.C*xI; 

            obj.bc.ct = zeros(obj.q,1);
            obj.bc.ctt = [0, 0, 0, 0, 0, 0;
                   0, 2*Dy1, 2*Dz1,0,0,0;
                   0, 0, 0, 0, 0, 0;
                   0, 2*Dy1-2*Dy2, 2*Dz1-2*Dz2, 0, 2*Dy2-2*Dy1, 2*Dz2-2*Dz1;
                   0, 0, 0, 0, 0, 0]*xII;
            
            
                       
        end
        function [xI_1, xI_2] = position(obj,y)
            xI_1 = y(:,2:3);
            xI_2 = y(:,5:6);
        end
        function [xII_1, xII_2] = velocity(obj,y,Dy)
            xII_1 = Dy(:,2:3);
            xII_2 = Dy(:,5:6);
        end
        function [alpha,beta] = angles(obj,y)
            alpha = wrapToPi(atan2(y(:,3),y(:,2))+pi/2)-pi/2;
            beta =  wrapToPi(atan2(y(:,3),y(:,5)-y(:,2))+pi/2)-pi/2;
        end
        
        function [Dalpha,Dbeta] = angleVelocities(obj,Dy)
            %% TODO
            Dalpha = Dy(:,1);
            Dbeta = Dy(:,4);
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
            c_f = matlabFunction(obj.bc.c,'vars',{[obj.xI;obj.xII]});
            Dc_f = matlabFunction(obj.bc.Dc,'vars',{[obj.xI;obj.xII]});
            
            T = length(t);
            
            c = zeros(T,length(obj.bc.c));
            Dc = zeros(T,length(obj.bc.Dc));
            for ii = 1:T        
                c(ii,:) = c_f([y(ii,:),Dy(ii,:)]');
                Dc(ii,:) = Dc_f([y(ii,:),Dy(ii,:)]');
            end
            
            results.c = c;
            results.Dc = Dc;
        end
    end
end

