%---Grid---%
M = 128;
T = 2*pi;
lambda = 0.25; %0.9,0.7071,0.25
dx = 2*pi/M;
x = 0: dx: 2*pi-dx;
dt = lambda*dx;

%---Initial Conditions---%
% IC = ones(1,length(x));
% IC(x<pi) = zeros(1,length(x(x<pi)));
k = 8;
IC = sin(k*x);
u = IC;

%---Lax-Wendroff---%
t = dt;
while t<T
    t = t+dt;
    u = u - lambda*([u(2:end) u(1)] - [u(end) u(1:end-1)])/2 + ...
        lambda^2*([u(2:end) u(1)] - 2*u + [u(end) u(1:end-1)])/2;
%     hold off
%     plot(x,u,'b')
%     hold on
%     plot(x,IC,'r')
%     pause(0.01)
end
hold on
plot(x,u,'m') %b,k,m
xlabel('x')
ylabel('u(x,t)')
set(gca,'FontSize',26)