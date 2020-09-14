function vis_results(cs,results,varargin)
        %% visualization parameters
    NN = length(results);
    cmap = colormap(lines(NN)); %colormap
    lwidth = 0.8;
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
    f = figure;
    for ii = 1:NN
        fig_2D{ii} = subplot(NN,1,ii,'Parent',f); hold on; grid on;
        xlabel("x [m]"); ylabel("y [m]");
        title("crankshaft in 2D (position of masses) " + results{ii}.task_info.name );
        daspect([1 1 1]);
        cb = colorbar;
        set(get(cb,"Label"),"string","velocity [m/s]");
    end
    % angles and angle speeds over time
    figure("Position",[680,325,563,653]); 
    fig_a1 = subplot(2,1,1); hold on; grid on;
    set(fig_a1,"outerposition",[0 0.65 1 0.35]);
    xlabel("t [s]"); ylabel("angles [rad]"); title("angles");
    fig_a2 = subplot(2,1,2); hold on; grid on;
    set(fig_a2,"outerposition",[0 0.0 1 0.65]);
    xlabel("t [s]"); ylabel("angle velocities [rad/s]"); title("angle velocities")
    legend('Location', 'southoutside' );
     
    % errors
     figure("Position",[680,325,563,653]);  
    for ii = 1:NN
        fig_e1{ii} = subplot(NN,2,2*ii-1); hold on; grid on;
    %     set(fig_e1, 'YScale', 'log');
        ylabel({results{ii}.task_info.name,"c(t)"});
        if(ii == 1)
            title("Position" );
        end
        if(ii == NN)
            xlabel("t [s]");
        end
%         legend('Location', 'bestoutside' );   

        fig_e2{ii} = subplot(NN,2,2*ii); hold on; grid on;
        xlabel("t [s]"); ylabel("dc(t)/dt");
        
        if(ii == 1)
            title("Velocity" );
        end
        if(ii == NN)
            xlabel("t [s]");
        end
%         legend('Location', 'bestoutside' );
    end
    try
        suptitle('Constraint Discrepancy')
    catch
    end
    for ii = 1:NN
       try
            % position of masses
            r1_ = r1(results{ii}.y);
            r2_ = r2(results{ii}.y);
            
            v1_ = v1([results{ii}.y,results{ii}.Dy]);
            v2_ = v2([results{ii}.y,results{ii}.Dy]);
            
            % errors of constraints
            plot(fig_e1{ii},cs.sym.tspan,results{ii}.c,'DisplayName',results{ii}.task_info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            plot(fig_e2{ii},cs.sym.tspan,results{ii}.Dc,'DisplayName',results{ii}.task_info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            
            % 2D positions
            t_step = 1;
            scatter(fig_2D{ii},r1_(1:t_step:end,1),r1_(1:t_step:end,2),[],v1_(1:t_step:end),'.','DisplayName',results{ii}.task_info.name);
            scatter(fig_2D{ii},r2_(1:t_step:end,1),r2_(1:t_step:end,2),[],v2_(1:t_step:end),'.','DisplayName',results{ii}.task_info.name);
            
            % angles and angle speeds over time           
            plot(fig_a1,cs.sym.tspan,results{ii}.y(:,1),'-','DisplayName',"\alpha - " + results{ii}.task_info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            plot(fig_a1,cs.sym.tspan,results{ii}.y(:,2),'--','DisplayName',"\beta - " + results{ii}.task_info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            plot(fig_a2,cs.sym.tspan,results{ii}.Dy(:,1),'-','DisplayName',"\alpha and d\alpha /dt - " + results{ii}.task_info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            plot(fig_a2,cs.sym.tspan,results{ii}.Dy(:,2),'--','DisplayName',"\beta and d\beta /dt - " + results{ii}.task_info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            
            %% animate 2D
            if ANIMATE
                animate_cs(cs.sym.tspan,r1_,r2_);
            end
       catch
          warning("something went wrong with result " + string(ii)); 
       end
    end
end