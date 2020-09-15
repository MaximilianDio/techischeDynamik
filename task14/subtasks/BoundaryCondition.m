%% Boundary condition: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Maximilian Dio - 21595892 - 
% Date:     15.09.2020
% Notes:    create derivatives, jacobian ct and ctt for implict boundary
% conditions
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef BoundaryCondition
    
    properties
        c
        Dc
        DDc
        
        C
        ct
        ctt
    end
    
    methods
        function obj = BoundaryCondition(c,y,Dy)
            obj.c = c;
            
            syms t_ real
            DDy = sym('DDy',size(y),'real');
            % taylor expansion
            y_ = y + t_*Dy + t_^2/2*DDy;
            
            % spans the space of fixed coordinates
            obj.C = jacobian(c,y);
            
            % differentiate bc
            c_ = subs(c,y,y_);
            Dc_ = diff(c_,t_);
            DDc_ = diff(Dc_,t_);
            
            obj.Dc = subs(Dc_,'t_',0);
            obj.DDc = subs(DDc_,'t_',0);
            
            obj.ct = simplify(obj.Dc - obj.C*Dy);
            obj.ctt = simplify(obj.DDc - obj.C*DDy);
            
        end
    end
end

