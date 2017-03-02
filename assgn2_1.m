%---Grid---%
M = 64;
dx = 2*pi/M;
x = 0:dx:2*pi-dx;
x = x';

T = 4*pi;
dt = 0.1;

% v = dx^2/(2*dt);
% a = dx*sqrt(1 - dt/2);
v = 0.1;
a = 1.155; %1.16 and 1.165


%%
%---Stability---%
figure(1)
alphas = v*2*(cos(dx*(1:M/2)) - 1)/(dx^2);
betas = a*sin(dx*(1:M/2))/dx;

% alpha = dt*v*2*(cos(dx*k) - 1)/(dx^2);
% beta = a*sin(dx*k)/dx;
% d1 = (1 - alpha*dt/2);
% b1 = -(1 + alpha*dt/2 + 1i*3*beta*dt/2);
% c1 = 1i*beta*dt/2;
% stab1 = max(abs(-b1 - sqrt(b1.^2 - 4*d1.*c1))./(2*d1),abs(-b1 + sqrt(b1.^2 - 4*d1.*c1))./(2*d1));
% stab1

[alphagrid , betagrid] = meshgrid([-100:0.05:0],[0:0.05:30]);
d = (1-alphagrid*dt/2);
b = -(1 + alphagrid*dt/2 + 1i*3*betagrid*dt/2);
c = 1i*betagrid*dt/2;
stab = max(abs(-b - sqrt(b.^2 - 4*d.*c))./(2*d),abs(-b + sqrt(b.^2 - 4*d.*c))./(2*d));

%Numerical stability region
pcolor(alphagrid,betagrid,abs(stab))
shading flat
colormap(jet)
colorbar
hold on
contour(alphagrid,betagrid,abs(stab),[1 1],'k','Linewidth',10)
% hold on
plot(alphas,betas,'r*','MarkerSize',10)
% plot(0.5*alphas,3*betas,'ro','MarkerSize',10)
xlabel('\alpha')
ylabel('\beta')
set(gca,'FontSize',26,'Linewidth',10)
hold off

%%
%---Stability for specific choice of a and v
for k = 1:M/2

[agrid, vgrid] = meshgrid([0:0.05:5],[0:0.05:5]);
alphagrid = vgrid*2*(cos(dx*(k)) - 1)/(dx^2);
betagrid = agrid*sin(dx*(k))/dx;
d = (1-alphagrid*dt/2);
b = -(1 + alphagrid*dt/2 + 1i*3*betagrid*dt/2);
c = 1i*betagrid*dt/2;
stab = max(abs(-b - sqrt(b.^2 - 4*d.*c))./(2*d),abs(-b + sqrt(b.^2 - 4*d.*c))./(2*d));

figure(2)
pcolor(agrid,vgrid,abs(stab))
shading flat
colormap(jet)
colorbar
hold on
contour(agrid,vgrid,abs(stab),[1 1],'k','Linewidth',10)
xlabel('a')
ylabel('v')
set(gca,'FontSize',26,'Linewidth',10)
plot(a,v,'k.','MarkerSize',10)
hold off
pause(0.01)

end

%%
%---Initial conditions---%
Initial = @(x) 1 - 2*heaviside(x-pi);
% Initial = @(x) sin(x).^2;
% k = 8;
% Initial = @(x) sin(k*x);
IC = Initial(x);
bc = Initial(0);
t = dt;

%Initialization of previous step
U_old = IC;
f_old = a*([U_old(2:end) ; U_old(1)] - [U_old(end) ; U_old(1:end-1)]);
g_old = [U_old(2:end) ; U_old(1)] - 2*U_old + [U_old(end) ; U_old(1:end-1)];
U_current = U_old + (dt/(2*dx))*f_old + (v*dt/dx^2)*g_old;

%Exact solution information
uhat_0 = fft(IC)/M;
kmode = [0:M/2-1 -M/2:-1]; kmode = kmode';
omega = a*1i*kmode - v*kmode.^2;

%---CNAB---$
y = v*dt/(2*dx^2);
A = gallery('tridiag',M,-y,1+2*y,-y);
A(1,end) = -y; A(end,1) = -y;
while t<T
    
    %Increment time step
    t = t+dt;
    
    %Run algorithm
    f_current = a*([U_current(2:end) ; U_current(1)] - [U_current(end) ; U_current(1:end-1)]);
    f_old = a*([U_old(2:end) ; U_old(1)] - [U_old(end) ; U_old(1:end-1)]);
    g_current = [U_current(2:end) ; U_current(1)] - 2*U_current + [U_current(end) ; U_current(1:end-1)];
    b = U_current + (dt/(4*dx))*(3*f_current - f_old) + y*g_current;
    
    U_new = A\b;
    U_old = U_current;
    U_current = U_new;
    
    %Exact solution
    uhat = uhat_0.*exp(omega*t);
    Exact = M*real(ifft(uhat));
    
    %Results after each time step
    figure(3)
    subplot(1,2,2)
    plot(x,U_current,'Linewidth',5)
    hold on
    plot(x,Exact,'r--','Linewidth',5)
    hold off
    xlabel('x')
    ylabel('u(x,t)')
    legend('CNAB','Exact')
    set(gca,'FontSize',26,'Linewidth',5)
    pause(0.01)
end