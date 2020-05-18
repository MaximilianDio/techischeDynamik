classdef Force < Drawable
    %FORCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vector
        endpoint
        attackpoint
        maxNorm
    end
    
    methods
        function obj = Force(vector,attackpoint)
            %FORCE vector and attacpoint
            obj.vector = vector;
            obj.attackpoint = attackpoint;
            
            obj.maxNorm = vecnorm(vector,2,1); 
        end
        
        function draw(obj,n)
             hold on
             
            % TODO calculate max 
             amplitude = norm(obj.vector(:,n));
             color = hsv2rgb([1 amplitude/obj.maxNorm 1]);
             vec = obj.vector(:,n)/amplitude;
    
            quiver3(obj.attackpoint(1,n),obj.attackpoint(2,n),obj.attackpoint(3,n),vec(1,n),vec(2,n),vec(3,n),0,...
                'MaxHeadSize',0.5,...
                'Color',color)
            hold off
        end

    end
end

