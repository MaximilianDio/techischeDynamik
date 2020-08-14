function animate_cs(t,r1,r2)
    margin = 0.1;

    figure; ax = axes; hold on;
    
    title("crankshaft");
    xmin = min([min(r1(:,1)),min(r2(:,1)),0]);
    xmax = max([max(r1(:,1)),max(r2(:,1)),0]);
    ymin = min([min(r1(:,2)),min(r2(:,2)),0]);
    ymax = max([max(r1(:,2)),max(r2(:,2)),0]);
    xlim([xmin-margin,xmax+margin]);
    ylim([ymin-margin,ymax+margin]);
    xlabel("x [m]"); ylabel("y [m]");
    daspect([1,1,1]);
    
    origin = zeros(size(r1));
    
    rod1 = Rod2D(ax,origin,r1,'r');
    rod2 = Rod2D(ax,r1,r2,'b');
    
    time_txt = text(xmin,ymin,"time: " + string(0.0) + " sec");
    
    
    for ii = 1:length(t)
        % update point
        rod1.draw(ii);
        rod2.draw(ii);
        
        time_txt.String = "time: " + string(t(ii) + " sec");
        
        pause(t(2)); 
    end
end
