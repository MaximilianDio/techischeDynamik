classdef BicycleWheel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        rOB_K; % [3,N] in m
        qZYX; %quaternions
        t_span;
        
        fig;
        ax;
        timeTxt
        
        cyl; % struct for cylinder coordinates
        com; % center of mass;
        
        cf; % body coordinate frame 
        cfK; % ground coordinate frame 
    end
    
    methods
        function obj = BicycleWheel(radius,thickness, position, q,t_span)
            
            
            obj.t_span = t_span;
            
            obj.qZYX = q;
            
            obj.rOB_K = position;
            
            obj.fig = figure('Name','BicycleWheel','Position',[2 42 958 954]);
            obj.ax = axes();
            
            % access parent folder
            addpath('..\');
            qInit = euler2quat([0,0,0]);
            
            %% create cylinder
            obj.cyl = Cylinder([0;0;0],qInit,radius,thickness,thickness);
            
            %% create location of COM
            obj.com = plot3(position(1,1),position(2,1),position(3,1),'r','LineWidth',1.5);
            
            %% create local coordinate frame coordinate 
            obj.cf = CoordinateFrame([0;0;0],qInit,0.2);
            
            obj.cfK = CoordinateFrame([0;0;0],euler2quat([0,0,0]),0.2);            
            %% time
            obj.timeTxt = text(1,1,1,["time" num2str(0) "seconds"],'FontSize', 14);
            
            
            %% define axis 
            title('Bicycle wheel gyroscope - dynamics','FontSize',20)
            
            axis([-1 1 -1 1 -1 1]);
            pbaspect([1 1 1]);
            xlabel('x in m','FontSize',12);
            ylabel('y in m','FontSize',12);
            zlabel('z in m','FontSize',12);
            
            obj.plot(1);
        end
        
        function plot(obj,n)
            % set time
            obj.timeTxt.String = ['time ' num2str(obj.t_span(n)) ' seconds'];
            
            % center of mass
            m = n;
            m = min(n,m);
            obj.com.XData = obj.rOB_K(1,n-m+1:n);
            obj.com.YData = obj.rOB_K(2,n-m+1:n);
            obj.com.ZData = obj.rOB_K(3,n-m+1:n);
            
            % wheel and coordinate frame
            obj.cf.transform(obj.rOB_K(:,n),obj.qZYX(:,n));
            obj.cyl.transform(obj.rOB_K(:,n),obj.qZYX(:,n));
        end
    end
end

