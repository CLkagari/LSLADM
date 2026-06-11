clc;tic;syms x t s;
a = 1;c = 2/5;

u0 = -8*c/3*cos(x/4)^2;u0x = diff(u0,x);u0xx = diff(u0x,x);u0xxx = diff(u0xx,x);

u1 = ilaplace(laplace(u0*u0xxx+u0*u0x+3*u0x*u0xx,t,s)/s^a,s,t);
u1x = diff(u1,x);u1xx = diff(u1x,x);u1xxx = diff(u1xx,x);

u2 = ilaplace(laplace((u1*u0xxx+u0*u1xxx)+(u1*u0x+u0*u1x)+3*(u1x*u0xx+u0x*u1xx),t,s)/s^a,s,t);
u2x = diff(u1,x);u2xx = diff(u1x,x);u2xxx = diff(u1xx,x);

u3 = ilaplace(laplace((u2*u0xxx+u0*u2xxx+u1*u1xxx)+(u2*u0x+u0*u2x+u1*u1x)+3*(u2x*u0xx+u0x*u2xx+u1x*u1xx),t,s)/s^a,s,t);

u = u0+u1+u2+u3;
ur = u0+1.0001*u1+0.9149*u2+0.0816*u3;
U = -8*c/3*cos(x/4-c*t/4)^2;
e = abs(u-U);ee = abs(ur-U);

r1 = symfun(ur,[t,x]);r2 = symfun(e,[t,x]);r3 = symfun(ee,[t,x]);
[t,x] = meshgrid(0:0.005:1,0:0.005:2*pi);
r4 = r1(t,x);r5 = r2(t,x);r6 = r3(t,x);
figure
mesh(t,x,r4);
xlabel('t');
ylabel('x')
zlabel('$\tilde{u}_3(x,t) , \alpha = 1$','Interpreter','latex');
figure
mesh(t,x,r5);
xlabel('t');
ylabel('x')
zlabel('$LADM error, \alpha = 1$','Interpreter','latex');
figure
mesh(t,x,r6);
xlabel('t');
ylabel('x')
zlabel('$LSLADM error, \alpha = 1$','Interpreter','latex');
toc;
