% Maximilian Dio 23.04.2020

classdef Rod2D < Drawable
    properties
        p1;
        p2; 
        rod; 
        link;
    end
    
    methods
        function obj = Rod2D(axes,p1,p2,color)
            obj.p1 = p1;
            obj.p2 = p2;
            
            % rod
            obj.rod = line(axes,[obj.p1(1,1),obj.p2(1,1)],[obj.p1(1,2),obj.p2(1,2)],'Color',color);
            % link
            obj.link = plot(axes,[obj.p1(1,1),obj.p2(1,1)],[obj.p1(1,2),obj.p2(1,2)],".",'MarkerSize',10,'Color',color);
        end
        
        function draw(obj,n)
            obj.rod.XData = [obj.p1(n,1),obj.p2(n,1)];
            obj.link.XData = [obj.p1(n,1),obj.p2(n,1)];
            
            obj.rod.YData = [obj.p1(n,2),obj.p2(n,2)];
            obj.link.YData = [obj.p1(n,2),obj.p2(n,2)];
        end
    end
end

