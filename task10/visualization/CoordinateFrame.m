classdef CoordinateFrame
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        X;
        Y; 
        Z;
        
        handle;
    end
    
    methods
        function obj = CoordinateFrame(origin,qInit,size)
            % vector of elementary coordinates and origin
        	X = [1 0 0 0]*size;
            Y = [0 1 0 0]*size; 
            Z = [0 0 1 0]*size;
            
            color = ['r','g','b'];
            
            % transform to initial pose
         	[obj.X,obj.Y,obj.Z] = quatTransformBFtoGF([X,0],[Y,0],[Z,0],origin,qInit);  
            
            obj.handle.line = gobjects(3,1);
            for ii = 1:3
                obj.handle.line(ii) = line([obj.X(4),obj.X(ii)],[obj.X(4),obj.Y(ii)],[obj.X(4),obj.Z(ii)],'Color',color(ii));
            end
        end
        function transform(obj,rOg,q)
            
            [Xl,Yl,Zl] = quatTransformBFtoGF(obj.X,obj.Y,obj.Z,rOg,q);
            for ii = 1:3
                obj.handle.line(ii).XData = [Xl(4),Xl(ii)];
                obj.handle.line(ii).YData = [Yl(4),Yl(ii)];
                obj.handle.line(ii).ZData = [Zl(4),Zl(ii)];
            end
        end
        
    end
end

