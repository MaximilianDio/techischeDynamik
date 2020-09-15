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
    
    %% seperate tree and 
    type_idx = zeros(NN,1);
    for ii = 1:NN
        type_idx(ii) = results{ii}.info.type == "tree" ;
    end
    NNtrees = sum(type_idx);
    NNredundant = NN-NNtrees;
    

    %% 2D representation of crankshaft movement
    f = figure("Position",[680,125,563,853]);  
    for ii = 1:NN
        fig_2D{ii} = subplot(NN,1,ii,'Parent',f); hold on; grid on;
        xlabel("x [m]"); ylabel("y [m]");
        title("crankshaft in 2D (position of masses) " + results{ii}.info.name );
        daspect([1 1 1]);
        cb = colorbar;
        set(get(cb,"Label"),"string","velocity [m/s]");
    end
    
    for ii = 1:NN
       try
            % position of masses
            [xI_1, xI_2] = cs{ii}.position(results{ii}.y);
            [xII_1, xII_2] = cs{ii}.velocity(results{ii}.y,results{ii}.Dy);
            
            % 2D positions
            t_step = 1;
            scatter(fig_2D{ii},xI_1(1:t_step:end,1),xI_1(1:t_step:end,2),[],(xII_1(1:t_step:end,1).^2+xII_1(1:t_step:end,2).^2).^(1/2),'.','DisplayName',results{ii}.info.name);
            scatter(fig_2D{ii},xI_2(1:t_step:end,1),xI_2(1:t_step:end,2),[],(xII_2(1:t_step:end,1).^2+xII_2(1:t_step:end,2).^2).^(1/2),'.','DisplayName',results{ii}.info.name);

            % animate 2D
            if ANIMATE
                animate_cs(results{ii}.t,xI_1,xI_2);
            end
       catch
          warning("something went wrong with result " + string(ii)); 
       end
    end
    
    %% angles and angle speeds over time
    figure("Position",[680,325,563,653]); 
    fig_a1 = subplot(2,1,1); hold on; grid on;
    yticks([-pi -pi/2 0 pi/2 pi])
    yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    set(fig_a1,"outerposition",[0 0.55 1 0.4]);
    xlabel("t [s]"); ylabel("angles [rad]"); title("angles");
    fig_a2 = subplot(2,1,2); hold on; grid on;
    set(fig_a2,"outerposition",[0 0.0 1 0.55]);
    xlabel("t [s]"); ylabel("angle velocities [rad/s]"); title("angle velocities")
    legend('Location', 'southoutside',"NumColumns",2 );
    for ii = 1:NN
        if results{ii}.info.type ~= "tree" 
            continue % skip if not tree structure
        end
        try
            [alpha, beta] = cs{ii}.angles(results{ii}.y);
            [Dalpha, Dbeta] = cs{ii}.angleVelocities(results{ii}.Dy);
            
            % angles and angle speeds over time           
            plot(fig_a1,results{ii}.t,alpha,'-','DisplayName',"\alpha - " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            plot(fig_a1,results{ii}.t,beta,'--','DisplayName',"\beta - " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            plot(fig_a2,results{ii}.t,Dalpha,'-','DisplayName',"\alpha and d\alpha /dt - " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            plot(fig_a2,results{ii}.t,Dbeta,'--','DisplayName',"\beta and d\beta /dt - " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);

       catch
          warning("something went wrong with result " + string(ii)); 
       end
    end
    
    %% errors tree
    figure("Position",[680,325,563,653]);  
    for ii = 1:NNtrees
        fig_e_tree{ii} = subplot(NNtrees,1,ii); hold on; grid on;
    %     set(fig_e1, 'YScale', 'log');
        ylabel({results{ii}.info.name,"c(t) and dc(t)/dt"});
        if(ii == NNtrees)
            xlabel("t [s]");
        end
        legend('Location','northeast');   

    end
    try
        suptitle('Constraint Discrepancy')
    catch
    end
    
    jj = 0;
    for ii = 1:NN
       if results{ii}.info.type == "tree" 
            jj = jj +1;
       else 
           continue;
       end
       try       
            % errors of constraints
            plot(fig_e_tree{jj},results{ii}.t,results{ii}.c,'DisplayName','c','Color',cmap(ii,:),'Linewidth',lwidth);
            plot(fig_e_tree{jj},results{ii}.t,results{ii}.Dc,'DisplayName','dc/dt','Color',cmap(ii,:),'LineStyle','-.','Linewidth',lwidth);
            
       catch
          warning("something went wrong with result " + string(ii)); 
       end
    end
    
        %% angles and angle speeds over time
    figure("Position",[680,325,563,453]); 
    fig_x{1} = subplot(2,2,1); hold on; grid on;
    xlabel("t [s]"); ylabel("x_{I,1} [m]"); title("positions mass 1");
    legend("Location","southoutside","NumColumns",2);
    fig_x{2} = subplot(2,2,2); hold on; grid on;
    xlabel("t [s]"); ylabel("x_{I,2} [m]"); title("positions mass 2");
    legend("Location","southoutside","NumColumns",2);
    fig_x{3} = subplot(2,2,3); hold on; grid on;
    xlabel("t [s]"); ylabel("x_{II,1} [m]"); title("velocities mass 1");
    legend("Location","southoutside","NumColumns",2);
    fig_x{4} = subplot(2,2,4); hold on; grid on;
    xlabel("t [s]"); ylabel("x_{II,2} [m]"); title("velocities mass 2");
    legend("Location","southoutside","NumColumns",2);
    for ii = 1:NN
        if results{ii}.info.type == "tree" 
            continue % skip if tree structure
        end
        try
            % position of masses
            [xI_1, xI_2] = cs{ii}.position(results{ii}.y);
            [xII_1, xII_2] = cs{ii}.velocity(results{ii}.y,results{ii}.Dy);
            
            % angles and angle speeds over time           
            plot(fig_x{1},results{ii}.t,xI_1(:,1),'-','DisplayName',"x_1 " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            plot(fig_x{1},results{ii}.t,xI_1(:,2),'--','DisplayName',"y_1 " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            
            plot(fig_x{2},results{ii}.t,xI_2(:,1),'-','DisplayName',"x_2 " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            plot(fig_x{2},results{ii}.t,xI_2(:,2),'--','DisplayName',"y_2 " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            
            plot(fig_x{3},results{ii}.t,xII_1(:,1),'-','DisplayName',"v_{x1} " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            plot(fig_x{3},results{ii}.t,xII_1(:,2),'--','DisplayName',"v_{y1} " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            
            plot(fig_x{4},results{ii}.t,xII_2(:,1),'-','DisplayName',"v_{x2} " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            plot(fig_x{4},results{ii}.t,xII_2(:,2),'--','DisplayName',"v_{y2} " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
            
       catch
          warning("something went wrong with result " + string(ii)); 
       end
    end
    
    %% errors redundant coordinates
    figure("Position",[680,325,563,453]);  
     jj = 0;
    for ii = 1:NN
        if results{ii}.info.type ~= "tree" 
            jj = jj +1;
        else 
            continue;
        end
        fig_e_redundant{1,jj} = subplot(NNredundant,2,2*jj-1); hold on; grid on;
        legend;

        ylabel({results{ii}.info.name,"c(t) and dc(t)/dt"});
        if(jj == 1)
            title("geometric constraints")
        end
        if(jj == NNredundant)
            xlabel("t [s]");
        end
        
        fig_e_redundant{2,jj} = subplot(NNredundant,2,2*jj); hold on; grid on;
        legend;
        if(jj == 1)
            title("velocity constraints")
        end
        if(jj == NNredundant)
            xlabel("t [s]");
        end
        

    end
    try
        suptitle('Constraint Discrepancy')
    catch
    end
    
    jj = 0;
    for ii = 1:NN
       if results{ii}.info.type ~= "tree" 
            jj = jj +1;
       else 
           continue;
       end
       try      
           lstyles = {"-","--","-","-.",":"};
            % errors of constraints
            for kk = 1:size(results{ii}.c,2)
                plot(fig_e_redundant{1,jj},results{ii}.t,results{ii}.c(:,kk),'DisplayName','c_' + string(kk),'Color',cmap(ii,:),'LineStyle',lstyles{kk},'Linewidth',lwidth);
                plot(fig_e_redundant{2,jj},results{ii}.t,results{ii}.Dc(:,kk),'DisplayName','Dc_' + string(kk),'Color',cmap(ii,:),'LineStyle',lstyles{kk},'Linewidth',lwidth);
            end

            
       catch
          warning("something went wrong with result " + string(ii)); 
       end
    end
    
    
end