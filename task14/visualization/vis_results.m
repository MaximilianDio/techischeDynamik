function vis_results(cs,results,varargin)

    close all
        %% visualization parameters
    NN = length(results);
    cmap = colormap(lines(NN)); %colormap    
    lwidth = 1;
    % flags
    ANIMATE = false;
    
    % plot flags
    POS_2D = true;
    ANGLES_OVER_TIME = true;
    ERRORS_TREE = true;
    REDUNDANT_OVER_TIME = true;
    ERRORS_REDUNDANT = true;
    COMPARE_ALL_2D = true;
    
    % save flag
    SAVE = false;
    SAVE_TXT = '';
    
    
        %% evaluate varargin 
    for ii =1:2:length(varargin)
        switch cell2mat(varargin(ii))
            case 'animate'
                ANIMATE = varargin{ii+1};
            case 'save'
                SAVE = varargin{ii+1};
            case 'text'
                SAVE_TXT = varargin{ii+1};
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
    if (POS_2D)
        f_pos2d = figure("Position",[680,125,563,853]);  
        for ii = 1:NN
            fig_2D{ii} = subplot(NN,1,ii,'Parent',f_pos2d); hold on; grid on;
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
    end
    %% angles and angle speeds over time
    if (ANGLES_OVER_TIME)
        f_angles = figure("Position",[680,325,563,653]); 

        fig_a1 = subplot(4,1,1); hold on; grid on;
        yticks([-pi -pi/2 0 pi/2 pi]); yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
        ylabel("angles [rad]"); title("angles");

        fig_a2 = subplot(4,1,2); hold on; grid on;
        ylabel("angle velocities [rad/s]"); title("angle velocities")

        fig_a3 = subplot(4,1,3); hold on; grid on;
        xlabel("t [s]"); ylabel("angle accelerations [rad/s^2]"); title("angle accelerations")

        % place legend in 4th subplot
        legend('Location', 'south',"NumColumns",2,"Position",[0.05 0.02 0.9 0.2]);

        for ii = 1:NN
            if results{ii}.info.type ~= "tree" 
                continue % skip if not tree structure
            end
            try
                [alpha, beta] = cs{ii}.angles(results{ii}.y);
                [Dalpha, Dbeta] = cs{ii}.angleVelocities(results{ii}.Dy);

                % plot over time
                % angles          
                plot(fig_a1,results{ii}.t,alpha,'-','DisplayName',"\alpha - " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                plot(fig_a1,results{ii}.t,beta,'--','DisplayName',"\beta - " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);

                % angle velocities
                plot(fig_a2,results{ii}.t,Dalpha,'-','DisplayName',"\alpha and d\alpha /dt - " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                plot(fig_a2,results{ii}.t,Dbeta,'--','DisplayName',"\beta and d\beta /dt - " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);

                % angle accelerations
                plot(fig_a3,results{ii}.t,results{ii}.DDy(:,1),'-','DisplayName',"\alpha, d\alpha/dt d^2\alpha/dt^2 - " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                plot(fig_a3,results{ii}.t,results{ii}.DDy(:,2),'--','DisplayName',"\beta, d\beta/dt d^2\beta/dt^2 - " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
           catch
              warning("something went wrong with result " + string(ii)); 
           end
        end
    end
    %% errors tree
    if(ERRORS_TREE)
        f_angles_error = figure("Position",[377,235,1096,653]);  
        for ii = 1:NNtrees
            % errors geometric
            fig_e{1,ii} = subplot(NNtrees,3,3*(ii-1)+1); hold on; grid on;
            ylabel({results{ii}.info.name,"c(t)"});
            if(ii == 1)
                title("c(t)");
            end
            if(ii == NNtrees)
                xlabel("t [s]");
            end
            % errors velocity
            fig_e{2,ii} = subplot(NNtrees,3,3*(ii-1)+2); hold on; grid on;
            ylabel("dc(t)/dt");
            if(ii == 1)
                title("dc(t)/dt");
            end
            if(ii == NNtrees)
                xlabel("t [s]");
            end
   
            % errors velocity
            fig_e{3,ii} = subplot(NNtrees,3,3*(ii-1)+3); hold on; grid on;
            ylabel("d^2c(t)/dt^2");
            if(ii == 1)
                title("d^2c(t)/dt^2");
            end
            if(ii == NNtrees)
                xlabel("t [s]");
            end

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
                plot(fig_e{1,jj},results{ii}.t,results{ii}.c,'DisplayName','c','Color',cmap(ii,:),'Linewidth',lwidth);
                plot(fig_e{2,jj},results{ii}.t,results{ii}.Dc,'DisplayName','dc/dt','Color',cmap(ii,:),'LineStyle','-.','Linewidth',lwidth);
                plot(fig_e{3,jj},results{ii}.t,results{ii}.DDc,'DisplayName','d^2c/dt^2','Color',cmap(ii,:),'LineStyle','-.','Linewidth',lwidth);
           catch
              warning("something went wrong with result " + string(ii)); 
           end
        end
    end
        %% redundant coordinates over time
    if(REDUNDANT_OVER_TIME)
        f_redundant = figure("Position",[897,314,947,637]);  
        for ii = 1:NN
            fig_x{1,1} = subplot(2,4,1); hold on; grid on;
            xlabel("t [s]"); ylabel("x_{I,1} [m]"); title("positions mass 1");
            fig_x{1,2} = subplot(2,4,2); hold on; grid on;
            xlabel("t [s]"); ylabel("x_{II,} [m/s]"); title("velocity mass 1");
            fig_x{1,3} = subplot(2,4,3); hold on; grid on;
            xlabel("t [s]"); ylabel("dx_{II,2}/dt [m/s^2]"); title("acceleration mass 1");
            legend("Position",[0.75,0.6,0.2,0.35]);
            
            fig_x{2,1} = subplot(2,4,4+1); hold on; grid on;
            xlabel("t [s]"); ylabel("x_{I,1} [m]"); title("positions mass 2");
            fig_x{2,2} = subplot(2,4,4+2); hold on; grid on;
            xlabel("t [s]"); ylabel("x_{II,} [m/s]"); title("velocity mass 2");
            fig_x{2,3} = subplot(2,4,4+3); hold on; grid on;
            xlabel("t [s]"); ylabel("dx_{II,2}/dt [m/s^2]"); title("acceleration mass 2");
            legend("Position",[0.75,0.1,0.2,0.35]);
        end

        for ii = 1:NN
            if results{ii}.info.type == "tree" 
                continue % skip if tree structure
            end
           try
                % position of masses
                [xI_1, xI_2] = cs{ii}.position(results{ii}.y);
                [xII_1, xII_2] = cs{ii}.velocity(results{ii}.y,results{ii}.Dy);
                [DxII_1, DxII_2] = cs{ii}.acceleration(results{ii}.y,results{ii}.Dy,results{ii}.DDy);
                  
                % angles and angle speeds over time           
                plot(fig_x{1,1},results{ii}.t,xI_1(:,1),'-','DisplayName',"x " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                plot(fig_x{1,1},results{ii}.t,xI_1(:,2),'--','DisplayName',"y " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                
                plot(fig_x{1,2},results{ii}.t,xII_1(:,1),'-','DisplayName',"x " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                plot(fig_x{1,2},results{ii}.t,xII_1(:,2),'--','DisplayName',"y " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                
              
                
                plot(fig_x{1,3},results{ii}.t,DxII_1(:,1),'-','DisplayName',"x_1 " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                plot(fig_x{1,3},results{ii}.t,DxII_1(:,2),'--','DisplayName',"y_1 " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);

                plot(fig_x{2,1},results{ii}.t,xI_2(:,1),'-','DisplayName',"x " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                plot(fig_x{2,1},results{ii}.t,xI_2(:,2),'--','DisplayName',"y " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                    
              
                
                plot(fig_x{2,2},results{ii}.t,xII_2(:,1),'-','DisplayName',"x " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                plot(fig_x{2,2},results{ii}.t,xII_2(:,2),'--','DisplayName',"y " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                
                plot(fig_x{2,3},results{ii}.t,DxII_2(:,1),'-','DisplayName',"x_2 " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                plot(fig_x{2,3},results{ii}.t,DxII_2(:,2),'--','DisplayName',"y_2 " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);

           catch
              warning("something went wrong with result " + string(ii)); 
           end
        end
    end
    %% errors redundant coordinates
    
    if(ERRORS_REDUNDANT)
        f_redundant_errors = figure("Position",[377,235,1096,653]); 
        jj = 0;
        for ii = 1:NN
            if results{ii}.info.type ~= "tree" 
                jj = jj +1;
            else 
                continue;
            end
            N = NNredundant;
            % errors geometric
            fig_e{1,jj} = subplot(N,3,3*(jj-1)+1); hold on; grid on;
            ylabel({results{ii}.info.name,"c(t)"});
            if(jj == 1)
                title("c(t)");
            end
            if(jj == N)
                xlabel("t [s]");
            end
            legend
            % errors velocity
            fig_e{2,jj} = subplot(N,3,3*(jj-1)+2); hold on; grid on;
            ylabel("dc(t)/dt");
            if(jj == 1)
                title("dc(t)/dt");
            end
            if(jj == N)
                xlabel("t [s]");
            end
            legend
            % errors acceleration
            fig_e{3,jj} = subplot(N,3,3*(jj-1)+3); hold on; grid on;
            ylabel("d^2c(t)/dt^2");
            if(jj == 1)
                title("d^2c(t)/dt^2");
            end
            if(jj == N)
                xlabel("t [s]");
            end
            legend;

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
               lstyles = {"-","--","-.",":"};
                % errors of constraints
                for kk = 1:size(results{ii}.c,2)
                    plot(fig_e{1,jj},results{ii}.t,results{ii}.c(:,kk),'DisplayName','c' + string(kk),'Color',cmap(ii,:),'LineStyle',lstyles{kk},'Linewidth',lwidth);
                    plot(fig_e{2,jj},results{ii}.t,results{ii}.Dc(:,kk),'DisplayName','dc' + string(kk) +'/dt ','Color',cmap(ii,:),'LineStyle',lstyles{kk},'Linewidth',lwidth);
                    plot(fig_e{3,jj},results{ii}.t,results{ii}.DDc(:,kk),'DisplayName','d^2c'+ string(kk) + '/dt^2 ','Color',cmap(ii,:),'LineStyle',lstyles{kk},'Linewidth',lwidth);
                end


           catch
              warning("something went wrong with result " + string(ii)); 
           end
        end
    end
   
    %% 2D representation of crankshaft movement
    if(COMPARE_ALL_2D)
        f_compare = figure("Position",[897,314,947,637]);  
        for ii = 1:NN
            fig_x{1,1} = subplot(2,4,1); hold on; grid on;
            xlabel("t [s]"); ylabel("x_{I,1} [m]"); title("positions mass 1");
            fig_x{1,2} = subplot(2,4,2); hold on; grid on;
            xlabel("t [s]"); ylabel("x_{II,1} [m/s]"); title("velocity mass 1");
            fig_x{1,3} = subplot(2,4,3); hold on; grid on;
            xlabel("t [s]"); ylabel("dx_{II,1}/dt [m/s^2]"); title("acceleration mass 1");
            legend("Position",[0.75,0.6,0.2,0.35]);
            
            fig_x{2,1} = subplot(2,4,4+1); hold on; grid on;
            xlabel("t [s]"); ylabel("x_{I,2} [m]"); title("positions mass 2");
            fig_x{2,2} = subplot(2,4,4+2); hold on; grid on;
            xlabel("t [s]"); ylabel("x_{II,2} [m/s]"); title("velocity mass 2");
            fig_x{2,3} = subplot(2,4,4+3); hold on; grid on;
            xlabel("t [s]"); ylabel("dx_{II,2}/dt [m/s^2]"); title("acceleration mass 2");
            legend("Position",[0.75,0.1,0.2,0.35]);
        end

        for ii = 1:NN
           try
                % position of masses
                [xI_1, xI_2] = cs{ii}.position(results{ii}.y);
                [xII_1, xII_2] = cs{ii}.velocity(results{ii}.y,results{ii}.Dy);
                [DxII_1, DxII_2] = cs{ii}.acceleration(results{ii}.y,results{ii}.Dy,results{ii}.DDy);
                  
                % angles and angle speeds over time           
                plot(fig_x{1,1},results{ii}.t,xI_1(:,1),'-','DisplayName',"x " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                plot(fig_x{1,1},results{ii}.t,xI_1(:,2),'--','DisplayName',"y " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                
                plot(fig_x{1,2},results{ii}.t,xII_1(:,1),'-','DisplayName',"x " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                plot(fig_x{1,2},results{ii}.t,xII_1(:,2),'--','DisplayName',"y " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                
              
                
                plot(fig_x{1,3},results{ii}.t,DxII_1(:,1),'-','DisplayName',"x_1 " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                plot(fig_x{1,3},results{ii}.t,DxII_1(:,2),'--','DisplayName',"y_1 " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);

                plot(fig_x{2,1},results{ii}.t,xI_2(:,1),'-','DisplayName',"x " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                plot(fig_x{2,1},results{ii}.t,xI_2(:,2),'--','DisplayName',"y " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                    
              
                
                plot(fig_x{2,2},results{ii}.t,xII_2(:,1),'-','DisplayName',"x " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                plot(fig_x{2,2},results{ii}.t,xII_2(:,2),'--','DisplayName',"y " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                
                plot(fig_x{2,3},results{ii}.t,DxII_2(:,1),'-','DisplayName',"x_2 " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);
                plot(fig_x{2,3},results{ii}.t,DxII_2(:,2),'--','DisplayName',"y_2 " + results{ii}.info.name,'Color',cmap(ii,:),'Linewidth',lwidth);

           catch
              warning("something went wrong with result " + string(ii)); 
           end
        end
    end
    %% save figures
    if (SAVE)
        saveas(f_pos2d,'figures/' + SAVE_TXT + '_pos2d.png');
        saveas(f_angles,'figures/' + SAVE_TXT + '_tree.png');
        saveas(f_angles_error,'figures/' + SAVE_TXT + '_tree_error.png');
        saveas(f_redundant,'figures/' + SAVE_TXT + '_redundant.png');
        saveas(f_redundant_errors,'figures/' + SAVE_TXT + '_redundant_error.png');
        saveas(f_compare,'figures/' + SAVE_TXT + '_compare.png');
    end
    
end