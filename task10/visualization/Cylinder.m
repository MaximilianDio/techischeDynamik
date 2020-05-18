classdef Cylinder
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        X;
        Y; 
        Z;
        
        origin;
        
        handle;
        star;
        stick;
    end
    
    methods
        function obj = Cylinder(origin,qInit,radius,height,thickness)
            
            %% create cylinder
            
            % different radii over height
            n = 5;
            r = [radius-thickness ; ones(n,1)*radius; radius-thickness];
            obj.origin = origin;
            % create cylinder
        	[X,Y,Z] = cylinder(r,50);
            Z = height*Z;
            
            % transform to initial pose
         	[obj.X,obj.Y,obj.Z] = quatTransformBFtoGF(X,Y,Z,origin,qInit); 
            
            % create surface plot
            obj.handle = surf(obj.X,obj.Y,obj.Z,'FaceColor','k','EdgeColor','none');
            hold on;
            obj.star = plot3(obj.X(5,1),obj.Y(5,1),obj.Z(5,1),'.','Color','r','MarkerSize',20);
            obj.stick = line([origin(1),origin(1)],[origin(2),origin(2)],[origin(3),origin(3)],'Color','k','LineWidth',2);
            
            
        end
        function transform(obj,rOg,q)
            
            [Xcyl,Ycyl,Zcyl] = quatTransformBFtoGF(obj.X,obj.Y,obj.Z,rOg,q);
            for ii = 1:3
                obj.handle.XData = Xcyl;
                obj.handle.YData = Ycyl;
                obj.handle.ZData = Zcyl;
            end
            
            obj.star.XData = Xcyl(5,1);
            obj.star.YData = Ycyl(5,1);
            obj.star.ZData = Zcyl(5,1);
            
            obj.stick.XData = [obj.origin(1),rOg(1)];
            obj.stick.YData = [obj.origin(2),rOg(2)];
            obj.stick.ZData = [obj.origin(3),rOg(3)];
        end
        
    end
end

