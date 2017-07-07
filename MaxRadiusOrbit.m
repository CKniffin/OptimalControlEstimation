function MaxRadiusOrbit
%Self-contained function that solves for the maximum radius orbit in a given
%time using Matlab's bvp4c and a shooting algorithm. The shooting algorithm
%may require a few minutes to solve.
% Written: C. Kniffin, 2016
    % Conversion Factors
    c_lb2kg = 0.453592;
    c_lb2n = 4.448222;
    c_km2au = 1.4959965e8;
    c_s2d = 60*60*24;
    
    % Given Parameters
    m0 = 10000*c_lb2kg; % Initial mass
    T = 0.85*c_lb2n*c_s2d^2/(c_km2au*1000);    %kgm/s^2 = c_s2d^2/(c_km2au*1000)kgAU/d^2
    dmdt = 12.9*c_lb2kg; % burn rate
    tf = 193; % 193 day orbit transfer
    r0 = 1; % initial orbit radius (AU)
    u0 = 0; % initial radial velocity
    mu = 1.3271244018e11*c_s2d^2/c_km2au^3;
    v0 = sqrt(mu/r0); % initial tangential velocity
    x0 = [r0 u0 v0];
    
    % BVP4C Solution
    solinit = bvpinit(linspace(0,tf,4),[r0 u0 v0 -1 -1 -1]);
    sol = bvp4c(@diffeq2,@bc2,solinit);
    t = linspace(0,tf,100);
    y = deval(sol,t);
    
    figure(1)
    subplot(3,1,1)
    plot(t,y(1,:))
    legend('r(t)')
    title('BVP4C Results')
    subplot(3,1,2)
    plot(t,y(2:3,:))
    legend('u(t)','v(t)')
    subplot(3,1,3)
    plot(t,atan2(y(5,:),y(6,:)).*180/pi+180)
    legend('\Phi (t)')
    xlabel('t (d)')
    
    figure(2)
    plot(t,y(4:6,:))
    legend('\lambda_r','\lambda_u','\lambda_v')
    title('BVP4C Results')
    xlabel('t (d)')
    
    % Shooting Algorithm Params
    tol = 1e-5;
    err = 1e6;
    kmax = 200;
    k = 1;
    alpha = 1;
    lambda0 = -[1;1;1];
%     lambda0 = [-1.877;-54;-117.7];
    lambda = lambda0;
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    
    %    
    while norm(err) > tol && k < kmax
        % Set ICs
%         in1 = reshape([0,0,0;0,0,0;0,0,1],1,9);
%         in2 = reshape([0,0,0;0,1,0;0,0,1],1,9);
        in1 = zeros(1,9);
        in2 = reshape(eye(3),1,9);
        inarg = [x0,lambda',in1,in2];
        % Solve diffeqs
        [tout,xout] = ode45(@diffeq,[0 tf],inarg,options);
        % Obtain error
        err = bc(xout(end,:)');
        % Obtain output sensitivity to input (du/dl0, dv/dl0, dlr/dl0)
        dxdl0 = reshape(xout(end,7:15),3,3);...
               %+[xout(end,7)/5,0,0;0,0,0;0,0,0];
               %+(norm(xout(end,7:15))/5)*eye(3);
               
        dpsidxf = [0 1 0;.5*sqrt(mu/(xout(end,1))^(3/2)) 0 1];%...
%                    -3*xout(end,6)*sqrt(mu)/(4*xout(end,1)^(5/2)),0,0];
               
%         dpsidl0 = dpsidxf*dxdl0;
%         dpsidl0 = dpsidl0 + norm(dpsidl0)*eye(3)/10;
        % Update ICs
        lambda = lambda - alpha*(pinv(dpsidxf*dxdl0)*err);
        % Increment Counter
        k=k+1;
        
        % Calculate Thrust Angle (Change range from -pi:pi to 0:2*pi)
        phi = atan2(xout(:,5),xout(:,6)).*180/pi;
        phi = phi + (phi<0)*360;
        disp(['Current error = ', num2str(norm(err))])
        figure(3)
        subplot(3,1,1)
        plot(tout,xout(:,1))    
        title('Shooter Results')
        legend('r (t)')
        subplot(3,1,2)
        plot(tout,xout(:,2:3))
        legend('u(t)','v(t)')
        subplot(3,1,3)
        plot(tout,phi)
        legend('\Phi (t)')
        xlabel('t (d)')
        figure(4)
        plot(tout,xout(:,4:6))
        legend('\lambda_r','\lambda_u','\lambda_v')
        title('Shooter Results')
        xlabel('t (d)')
        pause(.25)
    end
    
    disp(['Converged after ',num2str(k),' iterations'])
    
    figure(5)
    plot(xout(:,1),phi)
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
        ud = v^2/r - mu/r^2 + T*lu/sqrt(lu^2+lv^2)/(m0-dmdt*t);
        vd = -u*v/r + T*lv/sqrt(lu^2+lv^2)/(m0-dmdt*t);
        lrd = -lu*(-v^2/r^2+2*mu/r^3) - lv*(u*v/r^2);
        lud = -lr + lv*v/r;
        lvd = -lu*2*v/r + lv*u/r;
        
        dfdx = [0               ,     1,     0;...
               -v^2/r^2+2*mu/r^3,     0, 2*v/r;...
               u*v/r^2,            -v/r, -u/r];
        dfdl = [0,0,0;
                0,...
                -T/(m0-dmdt*t)*(lu^2/(lu^2 + lv^2)^(3/2) - 1/(lu^2 + lv^2)^(1/2)),...
                -T/(m0-dmdt*t)*((lu*lv)/(lu^2 + lv^2)^(3/2));...
                0,...
                -T/(m0-dmdt*t)*((lu*lv)/(lu^2 + lv^2)^(3/2)),...
                -T/(m0-dmdt*t)*(lv^2/(lu^2 + lv^2)^(3/2) - 1/(lu^2 + lv^2)^(1/2))];
                
       
        dhdx = [-lu*(2*v^2/r^3 - 6*mu/r^4) + 2*lv*u*v/r^3, -lv*v/r^2, 2*lu*v/r^2 - lv*u/r^2;...
                -lv*v/r^2, 0, lv/r; lu*2*v/r^2 - lv*u/r^2, lv/r, -2*lu/r];
                
        dhdl = [0, v^2/r^2 - 2*mu/r^3, -u*v/r^2;...
                -1, 0, v/r; 0, -2*v/r, u/r];

        
        z1d = dfdx*z1 + dfdl*z2;
        z2d = dhdx*z1 + dhdl*z2;

%         dfdm = [   1,     0, 0;...
%                    0, 2*v/r, 0;...
%                 -v/r,  -u/r, 0];
%                   
%         dfdn = [0,0,0;...
%                 -v^2/r^2+2*mu/r^3,...
%                 T/(m0-dmdt*t)*(lu^2/(lu^2 + lv^2)^(3/2) - 1/(lu^2 + lv^2)^(1/2)),...
%                 T/(m0-dmdt*t)*((lu*lv)/(lu^2 + lv^2)^(3/2));...
%                 u*v/r^2,...
%                 T/(m0-dmdt*t)*((lu*lv)/(lu^2 + lv^2)^(3/2)),...
%                 T/(m0-dmdt*t)*(lv^2/(lu^2 + lv^2)^(3/2) - 1/(lu^2 + lv^2)^(1/2))];
%             
%         dhdm = [-lv*v/r^2, 2*lu*v/r - lv*u/r^2, 0;...
%                 0, lv/r, 1;...
%                 lv/r, -lu*2/r, 0];
%         dhdn = [-lu*(2*v^2/r^3 - 6*mu/r^4) + 2*lv*u*v/r^3,v^2/r^2 - 2*mu/r^3, -u*v/r^2;...
%                 -lv*v/r^2, 0, v/r;...
%                 2*lu*v/r^2, -2*v/r, u/r];

%         dfdm = [0,2*v/r,0;...
%                 -v/r,-u/r,0;...
%                 -lv*v/r^2, 2*v*lu/r^2-lv*u/r^2, 0];
%             
%         dfdn = [-v^2/r^2+2*mu/r^3,...
%                 T/(m0-dmdt*t)*(lu^2/(lu^2 + lv^2)^(3/2) - 1/(lu^2 + lv^2)^(1/2)),...
%                 T/(m0-dmdt*t)*((lu*lv)/(lu^2 + lv^2)^(3/2));...
%                 u*v/r^2,...
%                 T/(m0-dmdt*t)*((lu*lv)/(lu^2 + lv^2)^(3/2)),...
%                 T/(m0-dmdt*t)*(lv^2/(lu^2 + lv^2)^(3/2) - 1/(lu^2 + lv^2)^(1/2));...
%                 -lu*(2*v^2/r^3 - 6*mu/r^4) + 2*lv*u*v/r^3,...
%                 v^2/r^2-2*mu/r^3, -u*v/r^2];
%                 
%         dhdm = [1, 0, 0;...
%                 0, lv/r, -1;...
%                 lv/r, -2*lu/r, 0];
%             
%         dhdn = [0,0,0;...
%                 -lv*v/r^2, 0, v/r;...
%                 2*lu*v/r^2-lv*u/r^2, -2*v/r, u/r];
%         
%         z1d = dfdm*z1 + dfdn*z2;
%         z2d = dhdm*z1 + dhdn*z2;
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
        rin = X0(1);
        uin = X0(2);
        vin = X0(3);
        rf = Xf(1);
        uf = Xf(2);
        vf = Xf(3);
        lrf = Xf(4);
        lvf = Xf(6);
        res = [rin-r0;uin-u0;vin-v0;...
            uf;vf-sqrt(mu/rf);lrf+1-lvf*sqrt(mu)/(2*(Xf(1))^(3/2))];
    end
    % Boundary Conditions
    function res = bc(Xf)
        rf = Xf(1);
        uf = Xf(2);
        vf = Xf(3);
        lrf = Xf(4);
        lvf = Xf(6);
        res = [uf;(vf-sqrt(mu/rf))];%;(lrf+1-lvf*sqrt(mu)/(2*(Xf(1))^(3/2)))];
    end
end