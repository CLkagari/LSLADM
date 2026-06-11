clc;tic;syms x t s;
a = 0.8;c = 2/5;

u0 = -8*c/3*cos(x/4)^2;u0x = diff(u0,x);u0xx = diff(u0x,x);u0xxx = diff(u0xx,x);

u1 = ilaplace(laplace(u0*u0xxx+u0*u0x+3*u0x*u0xx,t,s)/s^a,s,t);
u1x = diff(u1,x);u1xx = diff(u1x,x);u1xxx = diff(u1xx,x);

u2 = ilaplace(laplace((u1*u0xxx+u0*u1xxx)+(u1*u0x+u0*u1x)+3*(u1x*u0xx+u0x*u1xx),t,s)/s^a,s,t);
u2x = diff(u1,x);u2xx = diff(u1x,x);u2xxx = diff(u1xx,x);

u3 = ilaplace(laplace((u2*u0xxx+u0*u2xxx+u1*u1xxx)+(u2*u0x+u0*u2x+u1*u1x)+3*(u2x*u0xx+u0x*u2xx+u1x*u1xx),t,s)/s^a,s,t);

u = u0+u1+u2+u3;
ux = diff(u,x);uxx = diff(ux,x);uxxx =diff(uxx,x);
ur = u0+0.9991*u1+0.8926*u2+0.0991*u3;%a=0.8,solution
% ur = u0+0.9917*u1+0.8878*u2+0.0929*u3;%a=0.5,solution
% ur = u0+0.9740*u1+0.9321*u2+0.0356*u3;%a=0.2,solution
urx = diff(ur,x);urxx = diff(urx,x);urxxx =diff(urxx,x);

resu = abs((1*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^6)/(45*gamma(3/5)) - (1*t^(4/5)*gamma(13/5)/gamma(9/5)*sin(x/4)^6)/(45*gamma(3/5)) - (1*t^(4/5)*gamma(13/5)/gamma(9/5)*sin(x/4)^6)/(45*gamma(3/5)) + (1*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^6)/(45*gamma(3/5)) - (4*1*gamma(9/5)*cos(x/4)*sin(x/4)^3)/(15*gamma(4/5)) - (4*1*gamma(9/5)*cos(x/4)^3*sin(x/4))/(15*gamma(4/5)) - (2*1*t^(8/5)*gamma(17/5)/gamma(13/5)*cos(x/4)*sin(x/4)^7)/(315*gamma(2/5)) + (2*1*t^(8/5)*gamma(17/5)/gamma(13/5)*cos(x/4)^7*sin(x/4))/(315*gamma(2/5)) - (1*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^2*sin(x/4)^4)/(45*gamma(3/5)) + (1*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^4*sin(x/4)^2)/(45*gamma(3/5)) + (1*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^2*sin(x/4)^4)/(15*gamma(3/5)) + (1*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^4*sin(x/4)^2)/(9*gamma(3/5)) - (2*1*t^(8/5)*gamma(17/5)/gamma(13/5)*cos(x/4)^3*sin(x/4)^5)/(315*gamma(2/5)) + (2*1*t^(8/5)*gamma(17/5)/gamma(13/5)*cos(x/4)^5*sin(x/4)^3)/(315*gamma(2/5))-u*uxxx-u*ux-3*ux*uxx);%a=0.8,LADM residual
resur = abs((0.8926*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^6)/(45*gamma(3/5)) - (0.8926*t^(4/5)*gamma(13/5)/gamma(9/5)*sin(x/4)^6)/(45*gamma(3/5)) - (0.0991*t^(4/5)*gamma(13/5)/gamma(9/5)*sin(x/4)^6)/(45*gamma(3/5)) + (0.0991*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^6)/(45*gamma(3/5)) - (4*0.9991*gamma(9/5)*cos(x/4)*sin(x/4)^3)/(15*gamma(4/5)) - (4*0.9991*gamma(9/5)*cos(x/4)^3*sin(x/4))/(15*gamma(4/5)) - (2*0.0991*t^(8/5)*gamma(17/5)/gamma(13/5)*cos(x/4)*sin(x/4)^7)/(315*gamma(2/5)) + (2*0.0991*t^(8/5)*gamma(17/5)/gamma(13/5)*cos(x/4)^7*sin(x/4))/(315*gamma(2/5)) - (0.8926*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^2*sin(x/4)^4)/(45*gamma(3/5)) + (0.8926*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^4*sin(x/4)^2)/(45*gamma(3/5)) + (0.0991*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^2*sin(x/4)^4)/(15*gamma(3/5)) + (0.0991*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^4*sin(x/4)^2)/(9*gamma(3/5)) - (2*0.0991*t^(8/5)*gamma(17/5)/gamma(13/5)*cos(x/4)^3*sin(x/4)^5)/(315*gamma(2/5)) + (2*0.0991*t^(8/5)*gamma(17/5)/gamma(13/5)*cos(x/4)^5*sin(x/4)^3)/(315*gamma(2/5))-ur*urxxx-ur*urx-3*urx*urxx);%a=0.8,LSLADM residual
% resu = abs((8*1*t^(1/2)/gamma(3/2)*cos(x/4)^6)/375 + (8*1*t^(1/2)/gamma(3/2)*cos(x/4)^6)/375 - (8*1*t^(1/2)/gamma(3/2)*sin(x/4)^6)/375 - (8*1*t^(1/2)/gamma(3/2)*sin(x/4)^6)/375 - (8*1*t^(1/2)/gamma(3/2)*cos(x/4)^2*sin(x/4)^4)/375 + (8*1*t^(1/2)/gamma(3/2)*cos(x/4)^4*sin(x/4)^2)/375 + (8*1*t^(1/2)/gamma(3/2)*cos(x/4)^2*sin(x/4)^4)/125 + (8*1*t^(1/2)/gamma(3/2)*cos(x/4)^4*sin(x/4)^2)/75 - (32*1*gamma(3/2)*cos(x/4)*sin(x/4)^3)/(75*pi^(1/2)) - (32*1*gamma(3/2)*cos(x/4)^3*sin(x/4))/(75*pi^(1/2)) - (64*1*t*gamma(5/2)*cos(x/4)*sin(x/4)^7)/(5625*pi^(1/2)) + (64*1*t*gamma(5/2)*cos(x/4)^7*sin(x/4))/(5625*pi^(1/2)) - (64*1*t*gamma(5/2)*cos(x/4)^3*sin(x/4)^5)/(5625*pi^(1/2)) + (64*1*t*gamma(5/2)*cos(x/4)^5*sin(x/4)^3)/(5625*pi^(1/2))-u*uxxx-u*ux-3*ux*uxx);%a=0.5,LADM residual
% resur = abs((8*0.8878*t^(1/2)/gamma(3/2)*cos(x/4)^6)/375 + (8*0.0929*t^(1/2)/gamma(3/2)*cos(x/4)^6)/375 - (8*0.8878*t^(1/2)/gamma(3/2)*sin(x/4)^6)/375 - (8*0.0929*t^(1/2)/gamma(3/2)*sin(x/4)^6)/375 - (8*0.8878*t^(1/2)/gamma(3/2)*cos(x/4)^2*sin(x/4)^4)/375 + (8*0.8878*t^(1/2)/gamma(3/2)*cos(x/4)^4*sin(x/4)^2)/375 + (8*0.0929*t^(1/2)/gamma(3/2)*cos(x/4)^2*sin(x/4)^4)/125 + (8*0.0929*t^(1/2)/gamma(3/2)*cos(x/4)^4*sin(x/4)^2)/75 - (32*0.9917*gamma(3/2)*cos(x/4)*sin(x/4)^3)/(75*pi^(1/2)) - (32*0.9917*gamma(3/2)*cos(x/4)^3*sin(x/4))/(75*pi^(1/2)) - (64*0.0929*t*gamma(5/2)*cos(x/4)*sin(x/4)^7)/(5625*pi^(1/2)) + (64*0.0929*t*gamma(5/2)*cos(x/4)^7*sin(x/4))/(5625*pi^(1/2)) - (64*0.0929*t*gamma(5/2)*cos(x/4)^3*sin(x/4)^5)/(5625*pi^(1/2)) + (64*0.0929*t*gamma(5/2)*cos(x/4)^5*sin(x/4)^3)/(5625*pi^(1/2))-ur*urxxx-ur*urx-3*urx*urxx);%0.5,LSLADM residual
% resu = abs((4*1*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^6)/(75*gamma(2/5)) - (4*1*t^(1/5)*gamma(7/5)/gamma(6/5)*sin(x/4)^6)/(75*gamma(2/5)) - (4*1*t^(1/5)*gamma(7/5)/gamma(6/5)*sin(x/4)^6)/(75*gamma(2/5)) + (4*1*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^6)/(75*gamma(2/5)) - (16*1*gamma(6/5)*cos(x/4)*sin(x/4)^3)/(15*gamma(1/5)) - (16*1*gamma(6/5)*cos(x/4)^3*sin(x/4))/(15*gamma(1/5)) - (16*1*t^(2/5)*gamma(8/5)/gamma(7/5)*cos(x/4)*sin(x/4)^7)/(1125*gamma(3/5)) + (16*1*t^(2/5)*gamma(8/5)/gamma(7/5)*cos(x/4)^7*sin(x/4))/(1125*gamma(3/5)) - (4*1*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^2*sin(x/4)^4)/(75*gamma(2/5)) + (4*1*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^4*sin(x/4)^2)/(75*gamma(2/5)) + (4*1*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^2*sin(x/4)^4)/(25*gamma(2/5)) + (4*1*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^4*sin(x/4)^2)/(15*gamma(2/5)) - (16*1*t^(2/5)*gamma(8/5)/gamma(7/5)*cos(x/4)^3*sin(x/4)^5)/(1125*gamma(3/5)) + (16*1*t^(2/5)*gamma(8/5)/gamma(7/5)*cos(x/4)^5*sin(x/4)^3)/(1125*gamma(3/5))-u*uxxx-u*ux-3*ux*uxx);%a=0.2,LADM residual
% resur = abs((4*0.9321*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^6)/(75*gamma(2/5)) - (4*0.9321*t^(1/5)*gamma(7/5)/gamma(6/5)*sin(x/4)^6)/(75*gamma(2/5)) - (4*0.0356*t^(1/5)*gamma(7/5)/gamma(6/5)*sin(x/4)^6)/(75*gamma(2/5)) + (4*0.0356*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^6)/(75*gamma(2/5)) - (16*0.9740*gamma(6/5)*cos(x/4)*sin(x/4)^3)/(15*gamma(1/5)) - (16*0.9740*gamma(6/5)*cos(x/4)^3*sin(x/4))/(15*gamma(1/5)) - (16*0.0356*t^(2/5)*gamma(8/5)/gamma(7/5)*cos(x/4)*sin(x/4)^7)/(1125*gamma(3/5)) + (16*0.0356*t^(2/5)*gamma(8/5)/gamma(7/5)*cos(x/4)^7*sin(x/4))/(1125*gamma(3/5)) - (4*0.9321*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^2*sin(x/4)^4)/(75*gamma(2/5)) + (4*0.9321*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^4*sin(x/4)^2)/(75*gamma(2/5)) + (4*0.0356*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^2*sin(x/4)^4)/(25*gamma(2/5)) + (4*0.0356*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^4*sin(x/4)^2)/(15*gamma(2/5)) - (16*0.0356*t^(2/5)*gamma(8/5)/gamma(7/5)*cos(x/4)^3*sin(x/4)^5)/(1125*gamma(3/5)) + (16*0.0356*t^(2/5)*gamma(8/5)/gamma(7/5)*cos(x/4)^5*sin(x/4)^3)/(1125*gamma(3/5))-ur*urxxx-ur*urx-3*urx*urxx);%a=0.2,LSLADM residual

r1 = symfun(ur,[t,x]);r2 = symfun(resu,[t,x]);r3 = symfun(resur,[t,x]);
[t,x] = meshgrid(0:0.005:1,0:0.005:2*pi);
r4 = r1(t,x);r5 = r2(t,x);r6 = r3(t,x);
figure
mesh(t,x,r4);
xlabel('t');
ylabel('x')
zlabel('$\tilde{u}_3(x,t) , \alpha = 0.8$','Interpreter','latex');
figure
mesh(t,x,r5);
xlabel('t');
ylabel('x')
zlabel('$LADM residual, \alpha = 0.8$','Interpreter','latex');
figure
mesh(t,x,r6);
xlabel('t');
ylabel('x')
zlabel('$LSLADM residual, \alpha = 0.8$','Interpreter','latex');
toc;
