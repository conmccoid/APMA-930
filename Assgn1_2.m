%% Upstream Weighting with lambda = 1
%---Grid---%
lambda = 1;
M = 64;
T = 4*pi;
dx = 2*pi/M;
dt = lambda*dx;
x = 0: dx: 2*pi-dx;

%---Initial Conditions---%
t = dt;
IC = ones(1,length(x));
IC(x<pi) = zeros(1,length(x(x<pi)));
u = IC;

%---Upstream Weighting---%
figure(1)
subtightplot(2,2,1,[0.05,0.02])
while t<T
    t = t+dt;
    u = u - lambda*(u - [u(end) u(1:end-1)]);
    hold off
    plot(x,u,'b')
    hold on
    plot(x,IC,'r')
    pause(0.01)
end
title('Upstream weighting, \lambda = 1')
max(abs(u - IC))

%% Leap Frog with lambda = 1
%---Grid---%
lambda = 1;
M = 64;
T = 4*pi;
dx = 2*pi/M;
dt = lambda*dx;
x = 0: dx: 2*pi-dx;

%---Initial Conditions---%
t = dt+dt;
IC = ones(1,length(x));
IC(x<pi) = zeros(1,length(x(x<pi)));
w = IC;
v = w - lambda*(w - [w(end) w(1:end-1)]);
% v = w - 0.5*lambda*([w(2:end) w(1)] - [w(end) w(1:end-1)]);

%---Leap Frog and Centered Differences---%
subtightplot(2,2,2,[0.05,0.02])
while t<T
    t = t+dt;
    placeholder = w - lambda*([v(2:end) v(1)] - [v(end) v(1:end-1)]);
    w = v;
    v = placeholder;
    hold off
    plot(x,v,'b')
    hold on
    plot(x,IC,'r')
    pause(0.01)
end
title('Leap Frog, \lambda = 1')

%% Upstream Weighting with lambda = 0.5
%---Grid---%
lambda = 0.5;
M = 64;
T = 4*pi;
dx = 2*pi/M;
dt = lambda*dx;
x = 0: dx: 2*pi-dx;

%---Initial Conditions---%
t = dt;
IC = ones(1,length(x));
IC(x<pi) = zeros(1,length(x(x<pi)));
u = IC;

%---Upstream Weighting---%
subtightplot(2,2,3,[0.05,0.02])
while t<T
    t = t+dt;
    u = u - lambda*(u - [u(end) u(1:end-1)]);
    hold off
    plot(x,u,'b')
    hold on
    plot(x,IC,'r')
    pause(0.01)
end
title('Upstream Weighting, \lambda = 0.5')

%% Leap Frog with lambda = 0.5
%---Grid---%
lambda = 0.5;
M = 64;
T = 4*pi;
dx = 2*pi/M;
dt = lambda*dx;
x = 0: dx: 2*pi-dx;

%---Initial Conditions---%
t = dt+dt;
IC = ones(1,length(x));
IC(x<pi) = zeros(1,length(x(x<pi)));
w = IC;
% v = w;
v = w - 0.5*lambda*(w - [w(end) w(1:end-1)]);
% v = w - 0.5*lambda*([w(2:end) w(1)] - [w(end) w(1:end-1)]);

%---Leap Frog and Centered Differences---%
subtightplot(2,2,4,[0.05,0.02])
while t<T
    t = t+dt;
    placeholder = w - lambda*([v(2:end) v(1)] - [v(end) v(1:end-1)]);
    w = v;
    v = placeholder;
    hold off
    plot(x,v,'b')
    hold on
    plot(x,IC,'r')
    pause(0.01)
end
title('Leap Frog, \lambda = 0.5')