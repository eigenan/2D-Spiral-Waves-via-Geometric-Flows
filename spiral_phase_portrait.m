clc; close all; format longE;

% constants
V = 2;
D = 7;
omega = 1e-1; % list of omega values
alpha_star = zeros(size(omega)); % stores the alpha value on the trajectory with ell = 0
alpha_init = 1e-3:1e-2:1e-1; % list of initial values of alpha
ell_init = (-4e-2:1e-2:1e-1); % list of initial values of ell - omega/V
tspan = [1e+12 0]; % time interval for RK45

% stopping condition for RK45: halt whenever ell = 0
option1 = odeset('Events', @stop_condition);
option2 = odeset('RelTol',1e-8,'AbsTol',1e-10);
option = odeset(option1,option2);

% solve the ode for each omega value and store alpha_star as an array
for j=1:length(omega)
    for k=1:length(ell_init)
        for i=1:length(alpha_init)
            ell0 = omega(j)/V + ell_init(k);
            alpha0 = alpha_init(i);
            y0 = [ell0 alpha0]; % initial condition for RK45
            
            % solve the equation using RK45
            [t,y] = ode45(@(t,y) Fla(t,y,omega(j),D,V),tspan,y0,option);
            
            % store the alpha value with ell = 0 into alpha_star
            alpha_star(j) = y(end,end);

            if (abs(ell_init(k)) < 1e-4) && (i == 1)
                hold on
                % plot the solution on alpha-ell phase plane
                figure(j)
                plot(y(:,2),y(:,1),'-r','LineWidth',5)
                title('Phase Portrait with omega =',omega(j));
                xlabel('alpha');
                ylabel('ell');
            else
                hold on
                % plot the solution on alpha-ell phase plane
                figure(j)
                plot(y(:,2),y(:,1),'-b')
                title('Phase Portrait with omega =',omega(j));
                xlabel('alpha');
                ylabel('ell');
            % axis([0 y(end,2) min(y(:,1)) max(y(:,1))])
            end
        end
    end
end


% vector field
function dydt = Fla(t,y,omega,D,V)
    ell = y(1);
    alpha = y(2);
    dydt = zeros(2,1);
    dydt(1) = -(omega/D)*(alpha^2+ell^2)-2*(alpha^3)*ell-alpha*(ell^3)+(V/D)*(alpha^2+ell^2)^(3/2);
    dydt(2) = -alpha^4;
end

% stopping condition: ell = 0 or alpha > 10^4
function [value, isterminal, direction] = stop_condition(t,y)
    value      = [y(1) (1e+4 - y(2))]; % the values to reach zero
    isterminal = [1 1];   % stop integration
    direction  = [-1 -1]; % ell and 10^4 - alpha approach 0 from the positive side
end