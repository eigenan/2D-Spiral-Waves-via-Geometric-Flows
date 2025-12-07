close all; format longE;
file=load("data.mat");
font_size=18;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

alpha=file.alpha_star; % alpha_star array from simulation
omega = file.omega; % omega array used in the simulation
r=1./alpha;

% C0=2.338107410459763e+00; % smallest positive root of Bessel function
C0=1.018792971647471; % smallest positive root of Airy function
V=file.V;
D=file.D;

% figure(1)
% alpha_KS=omega./V + C0.*(2.^(1/3).*D.^(2/3).*V.^(-7/3)).*omega.^(5/3); % alpha KS
% polyfit(log(omega(1:10)),real(log(alpha(1:10) - alpha_KS(1:10))),1)
% plot(log(omega),real(log(alpha - alpha_KS)),'b.-')
% hold on
% title(['log(omega) vs log(alpha num - alpha KS) with V=',num2str(V),', D=',num2str(D)])
% xlabel('log(omega)')
% ylabel('log(alpha num - alpha KS)')
% % legend('results','Location','best')
% hold on

log_array=[real(log(r));real(log(omega-V./r))];
% for i=1:length(log_array(1,:))-1
%     if mod(i,2)~=0
%         if abs([log_array(1,i+1) log_array(2,i+1)]-[log_array(1,i) log_array(2,i)])<1
%             log_array(2,i)=NaN;
%             log_array(1,i)=NaN;
%         end
%     end
% end
figure(2)
plot(log_array(1,:),log_array(2,:),'ro','MarkerSize',6,'LineWidth',1)
hold on
% polyfit(log(r),real(log(omega-V./r)),1)
theo = - C0.*(2.^(1/3)).*(D.^(2/3)).*(V.^(1/3)).*(r.^(-5/3));
plot(real(log(r)),real(log(theo)),'-b','LineWidth',2)
hold on
ax=gca;
ax.FontSize = font_size-2;
xticks([1 3 5 7 9])
yticks([-14 -10 -6 -2])
title('$\log(R)$ vs $\log(\omega-VR^{-1})$','interpreter','latex','FontSize',font_size)
xlabel('$\log(R)$','interpreter','latex','FontSize',font_size)
ylabel('$\log(\omega-VR^{-1})$','interpreter','latex','FontSize',font_size)
legend('Numerical results','Asymptotic theory','interpreter','latex','FontSize',font_size)

% figure(4)
% plot(real(log(r)),real(log(abs(theo - (omega-V./r)))),'.r','MarkerSize',12)
% hold on
% polyfit(log(r),real(log(abs(theo-(omega-V./r)))),1)
% ax=gca;
% ax.FontSize = font_size-2;
% % xticks([1 3 5 7 9])
% % yticks([-14 -10 -6 -2])
% title('$\log(R)$ vs $\log($ numerical-theory $)$','interpreter','latex','FontSize',font_size)
% xlabel('$\log(R)$','interpreter','latex','FontSize',font_size)
% ylabel('$\log($difference$)$','interpreter','latex','FontSize',font_size)
% % legend('Numerical results','Asymptotic theory','interpreter','latex','FontSize',font_size)

figure(3)
plot(log(r(1:1:end)),omega(1:1:end),'ro','MarkerSize',6,'LineWidth',1)
hold on
omega_theo = V./r - C0.*(2.^(1/3)).*(D.^(2/3)).*(V.^(1/3)).*(r.^(-5/3));
plot(log(r(1:1:end)),omega_theo(1:1:end),'-b','LineWidth',2)
hold on
ax=gca;
ax.FontSize = font_size-2;
xticks([log(10) log(50) log(250) log(1000) log(5000)])
xticklabels({'10','50','250','1000','5000'})
yticks([0 0.05 0.1])
axis([0 10 0 0.1])
title('$R$ vs $\omega$','interpreter','latex','FontSize',font_size)
xlabel('$R$','interpreter','latex','FontSize',font_size)
ylabel('$\omega$','interpreter','latex','FontSize',font_size)
legend('Numerical results','Asymptotic theory','interpreter','latex','FontSize',font_size)