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

%     %% positions of masses
%     r1 = matlabFunction(cs.kin.r1,'vars',{[cs.dyn.yb]'});
%     r2 = matlabFunction(cs.kin.r2,'vars',{[cs.dyn.yb]'});
%     
%     v1 = matlabFunction(norm(cs.kin.v1),'vars',{[cs.dyn.yb;cs.dyn.Dyb]'});
%     v2 = matlabFunction(norm(cs.kin.v2),'vars',{[cs.dyn.yb;cs.dyn.Dyb]'});
    
    %%
    % 2D representation of crankshaft movement
    f = figure;
    for ii = 1:NN
        fig_2D{ii} = subplot(NN,1,ii,'Parent',f); hold on; grid on;
        xlabel("x [m]"); ylabel("y [m]");
        title("crankshaft in 2D (position of masses) " + results{ii}.info.name );
        daspect([1 1 1]);
        cb = colorbar;
        set(get(cb,"Label"),"string","velocity [m/s]");
    end
    % angles and angle speeds over time
    figure("Position",[680,325,563,653]); 
    fig_a1 = subplot(2,1,1); hold on; grid on;
    yticks([-pi -pi/2 0 pi/2 pi])
    yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    set(fig_a1,"outerposition",[0 0.65 1 0.35]);
    xlabel("t [s]"); ylabel("angles [rad]"); title("angles");
    fig_a2 = subplot(2,1,2); hold on; grid on;
    set(fig_a2,"outerposition",[0 0.0 1 0.65]);
    xlabel("t [s]"); ylabel("angle velocities [rad/s]"); title("angle velocities")
    legend('Location', 'southoutside' );
     
    % errors
     figure("Position",[680,325,563,653]);  
    for ii = 1:NN
        fig_e{ii} = subplot(NN,1,ii); hold on; grid on;
    %     set(fig_e1, 'YScale', 'log');
        ylabel({results{ii}.info.name,"c(t) and dc(t)/dt"});
        if(ii == NN)
            xlabel("t [s]");
        end
        legend('Location','northeast');   

    end
    try
        suptitle('Constraint Discrepancy')
    catch
    end
    for ii = 1:NN
       try
            % position of masses
            [xI_1, xI_2] = cs{ii}.position(results{ii}.y);
            [xII_1, xII_2] = cs{ii}.velocity(results{ii}.y,results{ii}.Dy);
            
            [alpha, beta] = cs{ii}.angles(results{ii}.y);
            [Dalpha, Dbeta] = cs{ii}.angleVelocities(results{ii}.Dy);
            
            % errors of constraints
            plot(fig_e{ii},results{ii}.t,results{ii}.c,'DisplayName','c','Color',cmap(ii,:),'Linewidth',lwidth);
            plot(fig_e{ii},results{ii}.t,results{ii}.Dc,'DisplayName','dc/dt','Color',cmap(ii,:),'LineStyle','-.','Linewidth',lwidth);
            
            % 2D positions
            t_step = 1;
            scatter(fig_2D{ii},xI_1(1:t_step:end,1),xI_1(1:t_step:end,2),[],(xII_1(1:t_step:end,1).^2+xII_1(1:t_step:end,2).^2).^(1/2),'.','DisplayName',results{ii}.info.name);
            scatter(fig_2D{ii},xI_2(1:t_step:end,1),xI_2(1:t_step:end,2),[],(xII_2(1:t_step:end,1).^2+xII_2(1:t_step:end,2).^2).^(1/2),'.','DisplayName',results{ii}.info.name);
            
            % angles and angle speeds over time           
            plot(fig_a1,results{ii}.t,alpha,'-','DisplayName',"\alpha - " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            plot(fig_a1,results{ii}.t,beta,'--','DisplayName',"\beta - " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            plot(fig_a2,results{ii}.t,Dalpha,'-','DisplayName',"\alpha and d\alpha /dt - " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            plot(fig_a2,results{ii}.t,Dbeta,'--','DisplayName',"\beta and d\beta /dt - " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            
            %% animate 2D
            if ANIMATE
                animate_cs(results{ii}.t,xI_1,xI_2);
            end
       catch
          warning("something went wrong with result " + string(ii)); 
       end
    end
end