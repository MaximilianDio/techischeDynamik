function vis_results(cs,results,varargin)
        %% visualization parameters
    % color
    col = ["b","r","g","m","k"];
    % flags
    ANIMATE = false;
    %% evaluate varargin 
    for ii =1:2:length(varargin)
        switch cell2mat(varargin(ii))
            case 'animate'
                ANIMATE = cell2mat(varargin(ii+1));
            otherwise
                error("parameter " +cell2mat(varargin(ii)) + " does not exsist");
        end
    end

    %% positions of masses
    r1 = matlabFunction(cs.kin.r1,'vars',{[cs.dyn.yb]'});
    r2 = matlabFunction(cs.kin.r2,'vars',{[cs.dyn.yb]'});
    
    v1 = matlabFunction(norm(cs.kin.v1),'vars',{[cs.dyn.yb;cs.dyn.Dyb]'});
    v2 = matlabFunction(norm(cs.kin.v2),'vars',{[cs.dyn.yb;cs.dyn.Dyb]'});
    
    %%
    % 2D representation of crankshaft movement
    figure; fig_2D = axes; hold on; grid on;
    xlabel("x [m]"); ylabel("y [m]");
    title("crankshaft in 2D (position of masses)");
    daspect([1 1 1]);
    cb = colorbar;
    set(get(cb,"Label"),"string","velocity [m/s]");
    
    % angles and angle speeds over time
    figure; 
    fig_a1 = subplot(2,1,1); hold on; grid on;
    xlabel("t [s]"); ylabel("angles [rad]"); title("angles");
    legend
    fig_a2 = subplot(2,1,2); hold on; grid on;
    xlabel("t [s]"); ylabel("angle velocities [rad/s]"); title("angle velocities")
    legend;
     % errors
    figure; 
    fig_e1 = subplot(2,1,1); hold on; grid on;
    xlabel("t [s]"); ylabel("c(t)");
    title("constraint discrepancy - position");
    legend;
    fig_e2 = subplot(2,1,2); hold on; grid on;
    xlabel("t [s]"); ylabel("dc(t)/dt");
    title("constraint discrepancy - velocity");
    legend;
    for ii = 1:length(results)
       try
            % position of masses
            r1_ = r1(results{ii}.y);
            r2_ = r2(results{ii}.y);
            
            v1_ = v1([results{ii}.y,results{ii}.Dy]);
            v2_ = v2([results{ii}.y,results{ii}.Dy]);
            
            % errors of constraints
            plot(fig_e1,cs.sym.tspan,results{ii}.c,'DisplayName',results{ii}.task_info.name,'Color',col(ii));
            plot(fig_e2,cs.sym.tspan,results{ii}.Dc,'DisplayName',results{ii}.task_info.name,'Color',col(ii));
            
            % 2D positions
            t_step = 1;
            scatter(fig_2D,r1_(1:t_step:end,1),r1_(1:t_step:end,2),[],v1_(1:t_step:end),'.','DisplayName',results{ii}.task_info.name);
            scatter(fig_2D,r2_(1:t_step:end,1),r2_(1:t_step:end,2),[],v2_(1:t_step:end),'.','DisplayName',results{ii}.task_info.name);
            
            % angles and angle speeds over time           
            plot(fig_a1,cs.sym.tspan,results{ii}.y(:,1),'-','DisplayName',"\alpha - " + results{ii}.task_info.name,'Color',col(ii));
            plot(fig_a1,cs.sym.tspan,results{ii}.y(:,2),'--','DisplayName',"\beta - " + results{ii}.task_info.name,'Color',col(ii));
            plot(fig_a2,cs.sym.tspan,results{ii}.Dy(:,1),'-','DisplayName',"d\alpha /dt - " + results{ii}.task_info.name,'Color',col(ii));
            plot(fig_a2,cs.sym.tspan,results{ii}.Dy(:,2),'--','DisplayName',"d\beta /dt - " + results{ii}.task_info.name,'Color',col(ii));
            
            %% animate 2D
            if ANIMATE
                animate_cs(cs.sym.tspan,r1_,r2_);
            end
       catch
          warning("something went wrong with result " + string(ii)); 
       end
    end
end