clc;tic;syms x y t s;
a = 1;k = 0.001;
u0 = 4/3*k*(sinh(x+y))^2;
u0x = diff(u0^2,x);u0xxx = diff(u0^2,x,3);u0xyy = diff(u0x,y,2);

u1 = -ilaplace(laplace(u0x+1/8*u0xxx+1/8*u0xyy,t,s)/s^a,s,t);
u1x = diff(2*u0*u1,x);u1xxx = diff(2*u0*u1,x,3);u1xyy = diff(u1x,y,2);

u2 = -ilaplace(laplace(u1x+1/8*u1xxx+1/8*u1xyy,t,s)/s^a,s,t);
u2x = diff(2*u0*u2+u1^2,x);u2xxx = diff(2*u0*u2+u1^2,x,3);u2xyy = diff(u2x,y,2);
 
u3 = -ilaplace(laplace(u2x+1/8*u2xxx+1/8*u2xyy,t,s)/s^a,s,t);
u = u0+u1+u2;
ur = u0+0.2897*u1-0.1997*u2;
U = 4/3*k*(sinh(x+y-k*t))^2;
e = abs(u-4/3*k*(sinh(x+y-k*t))^2);
ee = abs(ur-4/3*k*(sinh(x+y-k*t))^2);
t = 0.5;
r1 = symfun(ur,[t,x]);r2 = symfun(e,[t,x]);r3 = symfun(ee,[t,x]);
[x,y] = meshgrid(0:0.001:1,0:0.001:1);
r4 = r1(t,x);r5 = r2(t,x);r6 = r3(t,x);
figure
mesh(t,x,r4);
xlabel('t');
ylabel('x')
zlabel('$\tilde{u}_2(x,y,t) , \alpha = 1$','Interpreter','latex');
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