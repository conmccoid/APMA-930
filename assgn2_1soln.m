%---Grid---%
% M = 64;

Mn = 4:8; Mn = 2.^Mn;
Error = zeros(length(Mn),1);
for trial = 1:length(Mn)
    M = Mn(trial);

dx = 2*pi/M;
x = 0:dx:2*pi-dx;
x = x';

T = pi/2;
% dt = 0.1;
dt = dx;

v = 0.1;
a = 0.5;

% dtn = logspace(log10(4),log10(0.25),20); dtn = dtn*dx;
% 
% Error = zeros(length(dtn),1);
% 
% for trial = 1:length(dtn)
%     dt = dtn(trial);

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

[alphagrid , betagrid] = meshgrid(-100:0.05:0,0:0.05:30);
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
%---Initial conditions---%
% Initial = @(x) 1 - 2*heaviside(x-pi);
Initial = @(x) sin(x);
% k = 8;
% Initial = @(x) sin(k*x);
IC = Initial(x);
t = dt;

%Initialization of previous step
U_old = IC;
f_old = a*([U_old(2:end) ; U_old(1)] - [U_old(end) ; U_old(1:end-1)]);
g_old = [U_old(2:end) ; U_old(1)] - 2*U_old + [U_old(end) ; U_old(1:end-1)];
U_current = U_old + (dt/(2*dx))*f_old + (v*dt/dx^2)*g_old;

% %Exact solution information
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

    % Exact solution
%     Exact = exp(-v*t)*sin(x + a*t);
        
%     %Results after each time step
%     figure(3)
%     plot(x,U_current,'Linewidth',5)
%     hold on
%     plot(x,Exact,'r--','Linewidth',5)
%     hold off
%     xlabel('x')
%     ylabel('u(x,t)')
%     legend('CNAB','Exact')
%     set(gca,'FontSize',26,'Linewidth',5)
%     pause(0.01)
end

Error(trial) = norm(Exact - U_current)/norm(U_current);

end

% figure(2)
% loglog(1./dtn,Error,'Linewidth',2)
% % axis([0.9*dx 1.1*4*dx 1e-4 0.2])
% xlabel('1/\Delta t')
% ylabel('Relative error (2-norm)')
% set(gca,'FontSize',16,'Linewidth',2)

figure(2)
loglog(Mn,Error,'Linewidth',2)
% axis([0.9*dx 1.1*4*dx 1e-4 0.2])
xlabel('Number of grid points')
ylabel('Relative error (2-norm)')
set(gca,'FontSize',16,'Linewidth',2)