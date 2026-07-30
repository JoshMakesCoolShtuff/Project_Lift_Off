clear

th0 = 100;

zguess = [160;
          0.05;
          0.05;
          0];

ztrim = fsolve(@(z) trimfun(z,th0),zguess);

V0 = ztrim(1);
alpha0 = ztrim(2);
theta0 = ztrim(3);
el0 = ztrim(4);

x0 = [V0;
      alpha0;
      0;
      theta0;
      0;
      0];

u0 = [el0;
      th0];

fprintf('Trim values\n')
fprintf('V     = %.4f\n',V0)
fprintf('alpha = %.4f\n',alpha0)
fprintf('theta = %.4f\n',theta0)
fprintf('el    = %.4f\n',el0)

[A,B] = linmod_num(@CessnaEOM,x0,u0);

Along = A(1:4,1:4);
Blong = B(1:4,1);

pd = [-0.3+0.3i;
         -0.3-0.3i;
         -3;
         -4];

K = place(Along,Blong,pd)

Acl = Along - Blong*K;

dx0 = [5;
       0.03;
       0;
       0];

xi = x0;
xi(1:4) = x0(1:4) + dx0;

tspan = [0 100];

[t1,x1] = ode45(@(t,x) CessnaEOM(x,u0),tspan,xi);

[t2,x2] = ode45(@(t,x) cessnaClosed(t,x,x0,u0,K),tspan,xi);

el = zeros(length(t2),1);

for i = 1:length(t2)
    dx = x2(i,1:4)' - x0(1:4);
    el(i) = u0(1) - K*dx;
end

figure(1)
plot(t1,x1(:,1),'--',t2,x2(:,1),'LineWidth',1.5)
grid on
xlabel('t')
ylabel('V')
legend('Original','Feedback')
title('Velocity')

figure(2)
plot(t1,x1(:,2),'--',t2,x2(:,2),'LineWidth',1.5)
grid on
xlabel('t')
ylabel('\alpha')
legend('Original','Feedback')
title('Angle of Attack')

figure(3)
plot(t1,x1(:,3),'--',t2,x2(:,3),'LineWidth',1.5)
grid on
xlabel('t')
ylabel('q')
legend('Original','Feedback')
title('Pitch Rate')

figure(4)
plot(t1,x1(:,4),'--',t2,x2(:,4),'LineWidth',1.5)
grid on
xlabel('t')
ylabel('\theta')
legend('Original','Feedback')
title('Pitch Angle')

figure(5)
plot(t2,el,'LineWidth',1.5)
grid on
xlabel('t')
ylabel('e_l')
title('Elevator Deflection')

function F = trimfun(z,th0)

V = z(1);
alpha = z(2);
theta = z(3);
el = z(4);

x = [V;
     alpha;
     0;
     theta;
     0;
     0];

u = [el;
     th0];

xdot = CessnaEOM(x,u);

F = [xdot(1);
     xdot(2);
     xdot(3);
     theta - alpha];

end

function [A,B] = linmod_num(fun,x0,u0)

n = length(x0);
m = length(u0);

A = zeros(n,n);
B = zeros(n,m);

dx = 1e-5;
du = 1e-5;

for i = 1:n
    xp = x0;
    xm = x0;

    xp(i) = xp(i) + dx;
    xm(i) = xm(i) - dx;

    A(:,i) = (fun(xp,u0) - fun(xm,u0))/(2*dx);
end

for i = 1:m
    up = u0;
    um = u0;

    up(i) = up(i) + du;
    um(i) = um(i) - du;

    B(:,i) = (fun(x0,up) - fun(x0,um))/(2*du);
end

end

function xdot = cessnaClosed(t,x,x0,u0,K)

dx = x(1:4) - x0(1:4);

el = u0(1) - K*dx;
th = u0(2);

u = [el;
     th];

xdot = CessnaEOM(x,u);

end