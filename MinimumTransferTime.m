function MinimumTransferTime
%Self-contained function that solves for the minimum time to a given orbit
%(rf line 21). This code may take several minutes to converge. Needs some
%clean-up. To speed up run time the pause statement on line 103 may be
%removed.
% Written: C. Kniffin, 2016
close all;
    % Conversion Factors
    c_lb2kg = 0.453592;
    c_lb2n = 4.448222;
    c_km2au = 1.4959965e8;
    c_s2d = 60*60*24;
    
    % Given Parameters
    m0 = 10000*c_lb2kg; % initial spacecraft mass
    T = 0.85*c_lb2n*c_s2d^2/(c_km2au*1000);    %kgm/s^2 = c_s2d^2/(c_km2au*1000)kgAU/d^2
    dmdt = 12.9*c_lb2kg; % fuel burn rate
    tf = 195; % initial estimate for final time
    r0 = 1; % initial orbit (1AU)
    u0 = 0;
    mu = 1.3271244018e11*c_s2d^2/c_km2au^3;
    v0 = sqrt(mu/r0);
    x0 = [r0 u0 v0];
    rf = 1.537;
    uf = 0;
    vf = sqrt(mu/rf);
    
%     solinit = bvpinit(linspace(0,tf,4),[r0 u0 v0 -1 -1 -1]);
%     sol = bvp4c(@diffeq2,@bc2,solinit);
%     t = linspace(0,tf,100);
%     y = deval(sol,t);
%     
%     figure(1)
%     subplot(3,1,1)
%     plot(t,y(1,:))    
%     title('BVP4C Results')
%     subplot(3,1,2)
%     plot(t,y(2:3,:))
%     legend('u(t)','v(t)')
%     subplot(3,1,3)
%     plot(t,atan2(y(5,:),y(6,:)).*180/pi+180)
%     
%     figure(2)
%     plot(t,y(4:6,:))
%     legend('\lambda_r','\lambda_u','\lambda_v')
    
    % Shooting Algorithm Params
    tol = 1e-3; % Error tolerance
    kmax = 55; % maximum iterations
    k = 1; 
    alpha = .25; % gradient descent gain
    lambda0 = [1;1;1]; % initial guess for co-states
    % Upper limit and initial guess for final time
    tu = 200;
    tf = 195;
    % Lower Limit
    tl = 190;
    tprev = tl;
    lambda = lambda0;
    lambdagood = lambda;
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    itr = 1;
    
    while (abs(tprev - tf)) > .01
        disp(['Checking for tf = ',num2str(tf),'...'])
        % Reset error term
        err = 1e6*[1;1;1];
        while norm([err(1);err(2)*100;err(3)*100]) > tol && k < kmax
            % Set ICs
            in1 = zeros(1,9);
            in2 = reshape(eye(3),1,9);
            inarg = [x0,lambda',in1,in2];
            % Solve diffeqs
            [tout,xout] = ode45(@diffeq,[0 tf],inarg,options);
            % Obtain error
            err = bc(xout(end,:)');
            % Obtain output sensitivity to input (du/dl0, dv/dl0, dr/dl0)
            % and increase z1(1,1) by 20%
            grad = reshape(xout(end,7:15),3,3)...
                   + [xout(end,7)/5,0,0;0,0,0;0,0,0];%norm(xout(end,7:15))/10*eye(3);
            % Update ICs
            lambda = lambda - alpha*(grad\err);
            % Increment Counter
            k=k+1;
            disp(['Current error = ', num2str(norm([err(1);err(2)*100;err(3)*100]))])
            figure(3)
            subplot(3,1,1)
            plot(tout,xout(:,1))
            legend('r(t)')
            title('Shooter Results')
            subplot(3,1,2)
            plot(tout,xout(:,2:3))
            legend('u(t)','v(t)')
            subplot(3,1,3)
            plot(tout,atan2(xout(:,5),xout(:,6)).*180/pi+180)
            legend('\Phi(t)')
            xlabel('t (d)')
            figure(4)
            plot(tout,xout(:,4:6))
            legend('\lambda_r','\lambda_u','\lambda_v')
            xlabel('t (d)')
            title('Shooter Results')
            pause(.25)
%             tf = -T/dmdt*sqrt(xout(end,5)^2+xout(end,6)^2)+m0/dmdt;
        end
        if (norm(err) <= tol)
            disp(['Iteration ',num2str(itr),' converged after ',...
                num2str(k),' steps (tf = ',num2str(tf),')'])
            tprev = tf;
            tworking = tf;
            % Reset upper limit
            tu = tf;
            % Bisect 
            tf = (tl+tf)/2;
            % Store good lambda to use as initial guess 
            if (k > 1)
                lambdagood = lambda;
            end
        else
            disp(['Iteration ',num2str(itr),' did not converge (tf = ',num2str(tf),')'])
            tprev = tf;
            % Reset lower limit
            tl = tf;
            % Bisect
            tf = (tu+tf)/2;
            % Set lambda to last good lambda
            lambda = lambdagood;
        end
        k = 1;
        itr = itr + 1;
    end


%     err = 1e6*[1;1;1];        
%     k = 1;  
%     while norm([err(1);err(2)*100;err(3)*100]) > tol && k < kmax
%         % Set ICs
%         in1 = zeros(1,9);
%         in2 = reshape(eye(3),1,9);
%         inarg = [x0,lambda',in1,in2];
%         % Solve diffeqs
%         [tout,xout] = ode45(@diffeq,[0 tf],inarg,options);
%         % Obtain error
%         err = bc(xout(end,:)');
%         % Obtain output sensitivity to input (du/dl0, dv/dl0, dr/dl0)
%         % and increase z1(1,1) by 20%
%         grad = reshape(xout(end,7:15),3,3)...
%                + [xout(end,7)/5,0,0;0,0,0;0,0,0];%norm(xout(end,7:15))/10*eye(3);
%         % Update ICs
%         lambda = lambda - alpha*(grad\err);
%         % Increment Counter
%         k=k+1;
%         disp(['Current error = ', num2str(norm([err(1);err(2)*100;err(3)*100]))])
%         figure(3)
%         subplot(3,1,1)
%         plot(tout,xout(:,1))
%         legend('r(t)')
%         title('Shooter Results')
%         subplot(3,1,2)
%         plot(tout,xout(:,2:3))
%         legend('u(t)','v(t)')
%         subplot(3,1,3)
%         plot(tout,atan2(xout(:,5),xout(:,6)).*180/pi+180)
%         legend('\Phi(t)')
%         xlabel('t (d)')
%         figure(4)
%         plot(tout,xout(:,4:6))
%         legend('\lambda_r','\lambda_u','\lambda_v')
%         xlabel('t (d)')
%         title('Shooter Results')
%         pause(.25)
%         drhodtf = T*dmdt/(dmdt*tf-m0)^2*sqrt(xout(end,5)^2+xout(end,6)^2);
%         rho = -T/(m0-dmdt*tf)*sqrt(xout(end,5)^2+xout(end,6)^2);
%         tf = tf - alpha*rho/drhodtf
% %             tf = -T/dmdt*sqrt(xout(end,5)^2+xout(end,6)^2)+m0/dmdt;
%     end

    
    disp(['tf_min = ',num2str(tworking)])
    figure(5)
    plot(xout(:,1),atan2(xout(:,5),xout(:,6)).*180/pi+180)
    xlabel('r (AU)')
    ylabel('\Phi (^\circ)')
    title('Thrust Angle vs r')
    
    
    % Differential Equations
    function dxdt = diffeq(t,X)
        r = X(1);
        u = X(2);
        v=  X(3);
        lr = X(4);
        lu = X(5);
        lv = X(6);
        z1 = reshape(X(7:15),3,3);
        z2 = reshape(X(16:24),3,3);
        
        rd = u;
        ud = v^2/r - mu/r^2 - T*lu/sqrt(lu^2+lv^2)/(m0-dmdt*t);
        vd = -u*v/r - T*lv/sqrt(lu^2+lv^2)/(m0-dmdt*t);
        lrd = -lu*(-v^2/r^2+2*mu/r^3) - lv*(u*v/r^2);
        lud = -lr + lv*v/r;
        lvd = -lu*2*v/r + lv*u/r;
        
        dfdx = [0               ,     1,     0;...
               -v^2/r^2+2*mu/r^3,     0, 2*v/r;...
               u*v/r^2,            -v/r, -u/r];
        dfdl = [0,0,0;
                0,...
                T/(m0-dmdt*t)*(lu^2/(lu^2 + lv^2)^(3/2) - 1/(lu^2 + lv^2)^(1/2)),...
                T/(m0-dmdt*t)*((lu*lv)/(lu^2 + lv^2)^(3/2));...
                0,...
                T/(m0-dmdt*t)*((lu*lv)/(lu^2 + lv^2)^(3/2)),...
                T/(m0-dmdt*t)*(lv^2/(lu^2 + lv^2)^(3/2) - 1/(lu^2 + lv^2)^(1/2))];
                
       
        dhdx = [-lu*(2*v^2/r^3 - 6*mu/r^4) + 2*lv*u*v/r^3, -lv*v/r^2, 2*lu*v/r^2 - lv*u/r^2;...
                -lv*v/r^2, 0, lv/r; lu*2*v/r^2 - lv*u/r^2, lv/r, -2*lu/r];
                
        dhdl = [0, v^2/r^2 - 2*mu/r^3, -u*v/r^2;...
                -1, 0, v/r; 0, -2*v/r, u/r];
        
        z1d = dfdx*z1 + dfdl*z2;
        z2d = dhdx*z1 + dhdl*z2;

        dxdt = [rd;ud;vd;lrd;lud;lvd;reshape(z1d,9,1);reshape(z2d,9,1)];        
    end

    % Differential Equations (BVP4C)
    function dxdt = diffeq2(t,X)
        r = X(1);
        u = X(2);
        v=  X(3);
        lr = X(4);
        lu = X(5);
        lv = X(6);
        
        rd = u;
        ud = v^2/r - mu/r^2 - T*lu/sqrt(lu^2+lv^2)/(m0-dmdt*t);
        vd = -u*v/r - T*lv/sqrt(lu^2+lv^2)/(m0-dmdt*t);
        lrd = -lu*(-v^2/r^2+2*mu/r^3) - lv*(u*v/r^2);
        lud = -lr + lv*v/r;
        lvd = -lu*2*v/r + lv*u/r;

        dxdt = [rd;ud;vd;lrd;lud;lvd];        
    end

    % Boundary Conditions (BVP4C)
    function res = bc2(X0,Xf)
        res = [X0(1)-r0;X0(2)-u0;X0(3)-v0;Xf(1)-rf;Xf(2);Xf(3)-vf];
    end
    % Boundary Conditions
    function res = bc(Xf)
        res = [Xf(1)-rf;Xf(2);Xf(3)-vf];
    end
end