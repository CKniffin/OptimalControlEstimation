function SimultArrive
%Self-contained function that navigates two vehicles through an obstacle
%field such that they arrive at a goal location at the same time, while
%maintaining distance between each-other. Output is a gif that animates the
%path taken
% Written: C. Kniffin, 2016
    close all;
    global idx;
    % Vehicle Parameters
    V = 2; L = .25; maxphi = 40*pi/180;
    % Initial Conditions
    x10 = 3; y10 = 6; theta10 = -180*pi/180;
    x20 = 6; y20 = 4; theta20 = 10*pi/180;
    % Goal location
    xG = 7; yG = 15; 
    % Obstacles
    xObs = [10 9 7 6 5 4 3 2 1];
    yObs = [10 10 10 10 10 10 10 10 10];
    numObs = length(xObs); rObs = .5;    
    % ICs
    X0 = [x10 y10 theta10 x20 y20 theta20];
    X = X0; phi = [0 0];  
    % Time
    tf = 25; tspan = 0:.1:tf; tprev = tspan(1);
    % fmincon options
    options = odeset();
    optionsfmin = optimoptions(@fmincon,'Display','Off','Diagnostics','Off');
    ub = ones(2,1)*maxphi; lb = -ub;
    % Plot animation output
    filename = 'fmincon.gif';    
    % Gains
    kObs = 3; % Obstacle Gain
    
    for i = 2:length(tspan)
        tcurrent = tspan(i);
        idx = i;
        disp(['Current time = ',num2str(tcurrent)])
        phi(i,:) = fmincon(@costfun,phi(i-1,:),[],[],[],[],lb,ub,[],optionsfmin);
        X(i,:) = y(end,:);
        tprev = tcurrent;
        plotfun(i-1);
        if (y(end,1)-xG)^2+(y(end,2)-yG)^2<.5 && (y(end,4)-xG)^2+(y(end,5)-yG)^2<.5
            finaltime = tprev
            break
        end
    end
    
    figure
    plot(tspan(1:i),phi*180/pi)
    xlabel('Time')
    ylabel('\phi (deg)')
    legend('\phi_1','\phi_2')
    
    function J = costfun(phi)
    % Cost function J
        xvals = X(idx-1,:);
        pathxvals1 = X(:,1);
        pathxvals2 = X(:,4);
        pathyvals1 = X(:,2);
        pathyvals2 = X(:,5);
        % Simulate current time step
        [t,y]=ode45(@diffeq,[tprev tcurrent],xvals,options,phi);    
        
        x1 = y(end,1);
        y1 = y(end,2);
        x2 = y(end,4);
        y2 = y(end,5);
        % Distance to goal penalty
        distG1 = ((x1-xG)^2+(y1-yG)^2);
        distG2 = ((x2-xG)^2+(y2-yG)^2);
        % Ensure simultaneous arrival
        if distG1 < 3 && distG2 > 4
            JG1 = 1/distG1;
            JG2 = distG2;            
        elseif distG1 > 4 && distG2 < 3            
            JG1 = distG1;
            JG2 = 1/distG2;        
        else
            JG1 = distG1;
            JG2 = distG2;
        end
        JC = 0;
        % Loop over obstacles
        for n = 1:numObs
%             if sqrt((xObs(n)-x1)^2+(yObs(n)-y1)^2) < 1
                dist = sqrt((xObs(n)-x1)^2+(yObs(n)-y1)^2);
                JC = JC + exp(kObs/dist);
%             end
%             if sqrt((xObs(n)-x2)^2+(yObs(n)-y2)^2) < 1
                dist = sqrt((xObs(n)-x2)^2+(yObs(n)-y2)^2);
                JC = JC + exp(kObs/dist);
%             end
        end
        % Loop over discretized path
%         for n = 2:length(pathxvals1)-1
%             % Only add cost of distance is less than 2
%             if JG1 > 4
%                 dist = sqrt((pathxvals1(n)-x2)^2+(pathyvals1(n)-y2)^2);
%                 JC = JC + (30/dist);
%             end
%             if JG2 > 4
%                 dist = sqrt((pathxvals2(n)-x1)^2+(pathyvals2(n)-y1)^2);
%                 JC = JC + (30/dist);
%             end  
%         end
        if JG1 > 1 && JG2 > 1
            dist12 = sqrt((x1-x2)^2+(y1-y2)^2);
            JC = JC + 1/(dist12^3);
        end
        if JG1 < 2 && JG2 < 2
            JC = 0;
        end        
        J = 10*(JG1 + JG2) + JC;
    end

    function out = diffeq(t,x,phi)
        theta1 = x(3);
        theta2 = x(6);
        phi1 = phi(1);
        phi2 = phi(2);
%         phi1 = wrapTo2Pi(phi(1)); %Function not available in R2015a
%         phi2 = wrapTo2Pi(phi(2)); %Function not available in R2015a
        % Wrap phi to [-pi,pi]
        if phi1 > 2*pi
            phi1 = mod(phi1,2*pi);
            phi1 = phi1-2*pi;
        elseif phi1 < -2*pi
            phi1 = mod(phi1,2*pi);
            phi1 = phi1+2*pi;
        end
        if phi2 > 2*pi
            phi2 = mod(phi2,2*pi);
            phi2 = phi2-2*pi;
        elseif phi2 < -2*pi
            phi2 = mod(phi2,2*pi);
            phi2 = phi2+2*pi;
        end
        if phi1 > pi
            phi1 = phi1 - 2*pi;
        end
        if phi2 > pi
            phi2 = phi2 - 2*pi;
        end               
        
%         if phi1 > maxphi
%             phi1 = maxphi;
%         elseif phi1 < -maxphi
%             phi1 = -maxphi;
%         end
%         if phi2 > maxphi
%             phi2 = maxphi;
%         elseif phi2 < -maxphi
%             phi2 = -maxphi;
%         end
        
        xr1 = V*cos(phi1)*sin(theta1);
        xr2 = V*cos(phi2)*sin(theta2);
        yr1 = V*cos(phi1)*cos(theta1);
        yr2 = V*cos(phi2)*cos(theta2);
        theta1 = V/L*sin(phi1);
        theta2 = V/L*sin(phi2);
        out = [xr1 yr1 theta1 xr2 yr2 theta2]';
    end

    function plotfun(cntr)
        figure(1)
        h1 = plot(X(:,1),X(:,2),'b');
        hold on;
        h2 = plot(X(:,4),X(:,5),'r');
        for itr = 1:numObs
            h3 = plot(xObs(itr),yObs(itr),'x','color','black');
            circx = [linspace(-rObs,rObs,100)...
                     linspace(rObs,-rObs,100)];
            circy = sqrt(-circx(1:length(circx)/2).^2+rObs^2);
            circy(end+1:length(circx)) = -sqrt(-circx(1:length(circx)/2).^2+rObs^2);
            circx = circx + xObs(itr);
            circy = circy + yObs(itr);
            plot(circx,circy,'color','black')

        end
        circx = [linspace(-rObs,rObs,100)...
                 linspace(rObs,-rObs,100)];
        circy = sqrt(-circx(1:length(circx)/2).^2+rObs^2);
        circy(end+1:length(circx)) = -sqrt(-circx(1:length(circx)/2).^2+rObs^2);
        circx = circx + xG;
        circy = circy + yG;
        h4 = plot(xG,yG,'xg');
        plot(circx,circy,'g')
        legend([h1 h2 h3 h4],{'Vehicle 1','Vehicle 2','Obstacles','Goal'})
        axis equal
        ax = gcf;
        xlims = ax.CurrentAxes.XLim;
        ylims = ax.CurrentAxes.YLim;
        h = text(xlims(1)+(xlims(2)-xlims(1))/20,...
             ylims(2)-(ylims(2)-ylims(1))/20,...
             ['Current Time: ',num2str(tcurrent)]);
        if cntr ~= 0
            drawnow
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if cntr == 1;
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
              imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
            end
        end        
        delete(h);
    end
end