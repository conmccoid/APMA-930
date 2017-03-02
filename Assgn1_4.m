%---Grid---%
i = 3:10; i = 2.^i;
% i = [0.65, 0.50, 0.35, 0.1];
N = size(i);
U_error = zeros(N); H_error = zeros(N);
Error_est = zeros(N);
j = 0;
for M = i
    j = j+1;
    
% M = 64;
T = 2*pi;
lambda = 0.25;
G = 0.25;
dx = 2*pi/M;
dt = lambda*dx;
x = 0:dx:2*pi-dx;

%---Initial Conditions---%
% u0 = @(x) square(x);
% h0 = @(x) square(x - pi);
u0 = @(x) sin(x).^3;
h0 = @(x) cos(x).^3;
U0 = u0(x);
H0 = h0(x);
W = U0; K = H0;
U = U0 - lambda*(U0 - [U0(end) U0(1:end-1)]) - G*lambda*(H0 - [H0(end) H0(1:end-1)]);
H = H0 - lambda*(H0 - [H0(end) H0(1:end-1)]) - lambda*(U0 - [U0(end) U0(1:end-1)]);
f = @(x) (h0(x) + u0(x)/sqrt(G))/2; g = @(x) (-h0(x) + u0(x)/sqrt(G))/2;
% f = @(x,t) cos(x - t).*(cos(sqrt(G)*t) + sin(sqrt(G)*t)/sqrt(G));
% g = @(x,t) sqrt(G)*sin(x-t).*(cos(sqrt(G)*t)/sqrt(G) - sin(sqrt(G)*t));

%---Linearized Shallow Water Equations---%
t = dt;
while t<T
    t = t+dt;
    U_exact = sqrt(G)*(f(x - (1 + sqrt(G))*t) + g(x - (1 - sqrt(G))*t));
    H_exact = f(x - (1 + sqrt(G))*t) - g(x - (1 - sqrt(G))*t);
%     U_exact = g(x,t); H_exact = f(x,t);
    V = W - lambda*([U(2:end) U(1)] - [U(end) U(1:end-1)]) - ...
        G*lambda*([H(2:end) H(1)] - [H(end) H(1:end-1)]);
    J = K - lambda*([H(2:end) H(1)] - [H(end) H(1:end-1)]) - ...
        lambda*([U(2:end) U(1)] - [U(end) U(1:end-1)]);
    W = U; U = V;
    K = H; H = J;
    hold off
    plot(x,U_exact,'--b',x,H_exact,'--r')
    hold on
    plot(x,U,'b',x,H,'r')
    axis([0 2*pi -2.5 2.5])
%     plot(x,abs(U_exact - U),'b.',x,abs(H_exact - H),'r.')
    pause(0.01)
end

U_error(j) = norm(U_exact - U);
H_error(j) = norm(H_exact - H);
Error_est(j) = dt^2;

end

hold off
loglog(i,U_error,'.b',i,H_error,'.r',i,Error_est,'.k','MarkerSize',30)
xlabel('\lambda')
ylabel('Maximum error in u and h')
legend('Error in u','Error in h','dt^2')
set(gca,'FontSize',26)