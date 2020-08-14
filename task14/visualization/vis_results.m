function vis_results(cs,results)
    %% visualization parameters
    % color
    col = ["b","r","g","m","k"];
    
    %%
    r1 = matlabFunction(cs.kin.r1,'vars',{[cs.dyn.yb]'});
    r2 = matlabFunction(cs.kin.r2,'vars',{[cs.dyn.yb]'});
    
    % 2D representation of crankshaft movement
    figure; fig_2D = axes; hold on; grid on;
    xlabel("x [m]"); ylabel("y [m]");
    title("crankshaft in 2D (position of masses)");
    daspect([1 1 1]); 
    
    % angles and angle speeds over time
    figure; 
    fig_a1 = subplot(2,1,1); hold on; grid on;
    xlabel("t [s]"); ylabel("angles [rad]"); title("angles")
    fig_a2 = subplot(2,1,2); hold on; grid on;
    xlabel("t [s]"); ylabel("angle velocities [rad/s]"); title("angle velocities")
    
    for ii = 1:length(results)
       try
            r1_ = r1(results{ii}.y);
            r2_ = r2(results{ii}.y);
            
            
            
            % 2D
            plot(fig_2D,r1_(:,1),r1_(:,2),'-','Color',col(ii));
            plot(fig_2D,r2_(:,1),r2_(:,2),'-.','Color',col(ii));
            
            % angles and angle speeds over time           
            plot(fig_a1,cs.sym.tspan,results{ii}.y(:,1),'-','Color',col(ii));
            plot(fig_a1,cs.sym.tspan,results{ii}.y(:,2),'--','Color',col(ii));
            plot(fig_a2,cs.sym.tspan,results{ii}.Dy(:,1),'-','Color',col(ii));
            plot(fig_a2,cs.sym.tspan,results{ii}.Dy(:,2),'--','Color',col(ii));
            
            % errors of constraints
            figure; 
            subplot 211; hold on; grid on;
            plot(cs.sym.tspan,results{ii}.c);
            subplot 212; hold on; grid on;
            plot(cs.sym.tspan,results{ii}.Dc);
            
            %% animate 2D
            animate_cs(cs.sym.tspan,r1_,r2_);
       catch
          warning("something went wrong with result " + string(ii)); 
       end
    end
end