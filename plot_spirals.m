clc; close all; format longE;
font_size=24;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% constants
V = 1;
D = 1;
% omega = [2e-3 2e-2 2e-1 0.4]; % angular velocity
omega = [3.55e-2 2.1e-1];
% omega = 3.55e-2;
alpha_star = zeros(size(omega)); % stores the alpha value on the trajectory with ell = 0
alpha0 = 3.5e-3; % initial value of alpha
tspan = [1e+12 0]; % time interval for RK45

% stopping condition for RK45: halt whenever ell = 0
option1 = odeset('Events', @stop_condition);
option2 = odeset('RelTol',1e-8,'AbsTol',1e-10);
option = odeset(option1,option2);

% solve the ode for each omega value and store alpha_star as an array
for j=1:length(omega)

    ell0 = omega(j)/V + alpha0*D*omega(j)/(V^2); % initial value of ell
    h0 = 0;
    y0 = [h0 ell0 alpha0]; % initial condition for RK45
    
    % solve the equation using RK45
    [t,y] = ode45(@(t,y) Fhla(t,y,omega(j),D,V),tspan,y0,option);
    
    % alpha_star(j) = y(end,end); % store the alpha value with ell = 0 into alpha_star
    % phi = cumtrapz(flip(y(:,2))); % compute phi by integrating reversed ell
    % phi(floor(0.9.*length(phi)):end) = []; % trim data points near equilibrium
    % phi_periodic = mod(phi,2*pi); % make phi periodic
    % phi_periodic(abs(diff(phi_periodic))>2) = NaN; % remove vertical lines at jump points
    % r = flip(1./y(:,3));  % reverse the order of data because we integrated backward
    % r(floor(0.9.*length(r)):end) = []; % trim data points near equilibrium

    phi = flip(y(:,1)); % reverse phi
    % phi(floor(0.9.*length(phi)):end) = []; % trim data points near equilibrium
    phi_periodic = mod(phi,2*pi); % make phi periodic
    phi_periodic(abs(diff(phi_periodic))>2) = NaN; % remove vertical lines at jump points
    r = flip(1./y(:,3));  % reverse the order of data because we integrated backward
    % r(floor(0.9.*length(r)):end) = []; % trim data points near equilibrium

    % plot phi-r graph
    figure(j)
    plot(r,phi_periodic,'-b','LineWidth',2)
    title("$\phi$ - $r$ plot with core radius " + round(r(1)),'interpreter','latex','FontSize',font_size);
    xlabel('$r$','interpreter','latex','FontSize',font_size);
    ylabel('$\phi$','interpreter','latex','FontSize',font_size);
    ax=gca;
    ax.FontSize = font_size-2;
    % xticks([0 1000 2000 3000])
    yticks([0 pi 2*pi])
    yticklabels({'$0$','$\pi$','$2\pi$'})
    axis([0 200 0 2*pi])

    % plot the spiral
    figure(length(omega)+j)
    plot(r.*cos(phi),r.*sin(phi),'-r','LineWidth',2)
    hold on
    % core of radius r(1)
    angle_array = 0:pi/100:2*pi;
    plot(r(1).*cos(angle_array),r(1).*sin(angle_array),'-b','LineWidth',1)

    title("spiral with core radius " + round(r(1)),'interpreter','latex','FontSize',font_size);
    xlabel('$x$','interpreter','latex','FontSize',font_size);
    ylabel('$y$','interpreter','latex','FontSize',font_size);
    ax=gca;
    ax.FontSize = font_size-2;
    axis([-200 200 -200 200])
    axis square
    % xticks([-3000 -1500 0 1500 3000])
    % yticks([-3000 -1500 0 1500 3000])
    hold on
    
    % plot an inset on the first spiral plot
    if j == 1
        rectangle('Position',[-10 -60 25 25])
        ax2 = axes('Position',[0.58 0.66 0.25 0.25],'Box','on');
        plot(r.*cos(phi),r.*sin(phi),'-r','LineWidth',2)
        hold on
        plot(r(1).*cos(angle_array),r(1).*sin(angle_array),'-b','LineWidth',1)
        inset_pos = [-0.18 -0.04 -49.95 -49.8];
        axis(inset_pos)
        % ax2.FontSize = font_size-6;
        % axis([])
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
        axis square
    end
    get(get(gca,'title'),'position')

    
    



    % plot the spiral - the originl method
    % spiral_list = zeros(2,length(y(:,1)));
    % for k=1:length(spiral_list(1,:))
    %     spiral_list(:,k) = Spiral(y(k,1),y(k,3));
    % end
    % figure(2*length(omega)+j)
    % plot(spiral_list(1,:),spiral_list(2,:),'-r')
    % title('spiral with omega = ',omega(j));
    % xlabel('x');
    % ylabel('y');
    % hold on

    % plot the solution on alpha-ell phase plane
    % figure(2*length(omega)+j)
    % plot(y(:,3),y(:,2),'-b.')
    % title('Solution plot (alpha,ell) with $\omega=$',omega(j),'interpreter','latex');
    % xlabel('$\alpha$','interpreter','latex');
    % ylabel('$\ell$','interpreter','latex');

end

% spiral funciton in polar coordinates
% function spiral = Spiral(h,alpha)
%     spiral = zeros(2,1);
%     spiral(1) = (1/alpha)*cos(h);
%     spiral(2) = (1/alpha)*sin(h);
% end

% vector field
function dydt = Fhla(t,y,omega,D,V)
    ell = y(2);
    alpha = y(3);
    dydt = zeros(3,1);
    dydt(1) = ell*alpha^2;
    dydt(2) = -(ell^3)*alpha - 2*ell*(alpha^3) + (V/D)*(ell^2 + alpha^2)^(3/2) - sqrt(ell^2 + alpha^2)*sqrt(ell^2 ...
        + alpha^4)*(omega/D);
    dydt(3) = -alpha^4;
end

% stopping condition: ell = 0 or alpha > 10^4
function [value, isterminal, direction] = stop_condition(t,y)
    value      = [y(2) (1e+4 - y(3))]; % the values to reach zero
    isterminal = [1 1];   % stop integration
    direction  = [-1 -1]; % ell and 10^4 - alpha approach 0 from the positive side
end