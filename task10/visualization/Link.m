% Maximilian Dio 23.04.2020

classdef Link < Drawable
    properties ( Access = private)
       validColors = {'k','b','r','g','y'} 
    end
    properties
        m; % mass
        origin;
        r; % location vector r=r(y), y generalized coordinate
        color;
    end
    
    methods
        function obj = Link(m,r,origin,color)
            obj.color = color; %obj.validColors{randi([1,length(obj.validColors)])};
            obj.m = m;
            obj.r = r;
            obj.origin = origin;
        end
        
        function draw(obj,n)
            hold on
           % plots dot andl line to origin in current axis where n is the time-STEP 
           try
                line([obj.origin(1,n),obj.r(1,n)],[obj.origin(2,n),obj.r(2,n)],[obj.origin(3,n),obj.r(3,n)],'Color',obj.color);
           catch
                line([obj.origin(1),obj.r(1,n)],[obj.origin(2),obj.r(2,n)],[obj.origin(3),obj.r(3,n)],'Color',obj.color);
           end
            % link
            plot3(obj.r(1,n),obj.r(2,n),obj.r(3,n),".",'MarkerSize',10,'Color',obj.color);
            % travelled 
            plot3(obj.r(1,1:n),obj.r(2,1:n),obj.r(3,1:n),'-.','MarkerSize',2,'Color',obj.color);
            hold off
        end
    end
end

