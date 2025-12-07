R=3;
D=1;
V=1;

dx=.1;
L=100;
r=[R:dx:R+L]';
N=length(r);
psi=[0:.01:2*pi]';
%%%%%%%%%%%%%%%
% second order Laplacian with Dirichlet
%%%%%%%%%%%%%%%%%%%
e=ones(N,1);
lap=spdiags([e -2*e e],-1:1,N,N);
%%%%%%%%%%%%%%%%%%%%%
% Neumann
%%%%%%%%%%%%%%%%%%%%%
lap(1,1)=-1;
lap(N,N)=-1;
%
lap=lap/dx^2;
D2=lap/dx^2;

rdiag=spdiags(r,0,N,N);


D2r=spdiags(r.^(-2),0,N,N)*D2;

%
% %%%%%%%%%%%%%%%%%%%
% % centered FD second order
% %%%%%%%%%%%%%%%%%%%
% D1=spdiags([-1/2*e 0*e 1/2*e],-1:1,N,N);
% %%%%%%%%%%%%%%%%%%%%%
% % Neumann
% %%%%%%%%%%%%%%%%%%%%%
% D1(1,1)=-1/2;
% D1(N,N)=1/2;
% %
% D1=D1/dx;


%%%%%%%%%%%%%%%%%%%
% upwind
% backwardFD second order 2	−3/2	2	−1/2
%%%%%%%%%%%%%%%%%%%

D1=-spdiags([-1/2*e 2*e -3/2*e ],-2:0,N,N);
D1(1,1)=0; %Neumann
D1(2,1)=-3/2;%Neumann
D1=D1/dx;





%%%%%%%%%%%%%%%%%%%%5


phi=zeros(N,1);
phi=atan(r);
% phi=.1*sech(r-R);
% phi=-.01*r;
%%%%%%%%%%%%%%%%
dt=.003;
T=1000;
t=0;
t_plot=3;
h=figure(1);
A=speye(N,N)-dt*D2r;

while t<T
    t=t+dt;
    M=(1+r.^2.*(D1*phi).^2);
    F=-(D2*phi).*((D1*phi).^2)./M;
    F=F+(D1*phi).*(ones(N,1)+M)./(r.*M);
    F=F-V*M.^(1/2)./r;

%     F=0*phi;
    phin=A\(phi+dt*F);
    phi=phin;


    if mod(t,t_plot)<= dt % plotting, where most of the time is spent
        set(0,'CurrentFigure',h)
%         plot(r,mod(phi,2*pi),'b-')
        plot(r.*cos(phi),r.*sin(phi),'b-');
        hold on
        plot(R*cos(psi),R*sin(psi),'k-')
        plot((R+L)*cos(psi),(R+L)*sin(psi),'k-')
        hold off
        daspect([1 1 1])
        axis([-100 100 -100 100])
%         axis([r(1) r(N) -10 10])
        title(['t=' num2str(t) ',   max-min=' num2str(max(phi)-min(phi))])
        drawnow
    end
end



