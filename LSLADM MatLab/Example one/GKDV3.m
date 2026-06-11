syms x t s e1 e2 e3;
a = 0.8;c = 2/5;

u0 = -8*c/3*cos(x/4)^2;u0x = diff(u0,x);u0xx = diff(u0x,x);u0xxx = diff(u0xx,x);

u1 = ilaplace(laplace(u0*u0xxx+u0*u0x+3*u0x*u0xx,t,s)/s^a,s,t);
u1x = diff(u1,x);u1xx = diff(u1x,x);u1xxx = diff(u1xx,x);

u2 = ilaplace(laplace((u1*u0xxx+u0*u1xxx)+(u1*u0x+u0*u1x)+3*(u1x*u0xx+u0x*u1xx),t,s)/s^a,s,t);
u2x = diff(u1,x);u2xx = diff(u1x,x);u2xxx = diff(u1xx,x);

u3 = ilaplace(laplace((u2*u0xxx+u0*u2xxx+u1*u1xxx)+(u2*u0x+u0*u2x+u1*u1x)+3*(u2x*u0xx+u0x*u2xx+u1x*u1xx),t,s)/s^a,s,t);

u = u0+e1*u1+e2*u2+e3*u3;ux = diff(u,x);uxx = diff(ux,x);uxxx =diff(uxx,x);

%u11 = (diff(u,t)+diff(u^2,x)+1/8*diff(u^2,x,3)+1/8*diff(ux,y,2))^2;
ur08 = ((e2*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^6)/(45*gamma(3/5)) - (e2*t^(4/5)*gamma(13/5)/gamma(9/5)*sin(x/4)^6)/(45*gamma(3/5)) - (e3*t^(4/5)*gamma(13/5)/gamma(9/5)*sin(x/4)^6)/(45*gamma(3/5)) + (e3*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^6)/(45*gamma(3/5)) - (4*e1*gamma(9/5)*cos(x/4)*sin(x/4)^3)/(15*gamma(4/5)) - (4*e1*gamma(9/5)*cos(x/4)^3*sin(x/4))/(15*gamma(4/5)) - (2*e3*t^(8/5)*gamma(17/5)/gamma(13/5)*cos(x/4)*sin(x/4)^7)/(315*gamma(2/5)) + (2*e3*t^(8/5)*gamma(17/5)/gamma(13/5)*cos(x/4)^7*sin(x/4))/(315*gamma(2/5)) - (e2*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^2*sin(x/4)^4)/(45*gamma(3/5)) + (e2*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^4*sin(x/4)^2)/(45*gamma(3/5)) + (e3*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^2*sin(x/4)^4)/(15*gamma(3/5)) + (e3*t^(4/5)*gamma(13/5)/gamma(9/5)*cos(x/4)^4*sin(x/4)^2)/(9*gamma(3/5)) - (2*e3*t^(8/5)*gamma(17/5)/gamma(13/5)*cos(x/4)^3*sin(x/4)^5)/(315*gamma(2/5)) + (2*e3*t^(8/5)*gamma(17/5)/gamma(13/5)*cos(x/4)^5*sin(x/4)^3)/(315*gamma(2/5))-u*uxxx-u*ux-3*ux*uxx)^2;
%ur05 = ((8*e2*t^(1/2)/gamma(3/2)*cos(x/4)^6)/375 + (8*e3*t^(1/2)/gamma(3/2)*cos(x/4)^6)/375 - (8*e2*t^(1/2)/gamma(3/2)*sin(x/4)^6)/375 - (8*e3*t^(1/2)/gamma(3/2)*sin(x/4)^6)/375 - (8*e2*t^(1/2)/gamma(3/2)*cos(x/4)^2*sin(x/4)^4)/375 + (8*e2*t^(1/2)/gamma(3/2)*cos(x/4)^4*sin(x/4)^2)/375 + (8*e3*t^(1/2)/gamma(3/2)*cos(x/4)^2*sin(x/4)^4)/125 + (8*e3*t^(1/2)/gamma(3/2)*cos(x/4)^4*sin(x/4)^2)/75 - (32*e1*gamma(3/2)*cos(x/4)*sin(x/4)^3)/(75*pi^(1/2)) - (32*e1*gamma(3/2)*cos(x/4)^3*sin(x/4))/(75*pi^(1/2)) - (64*e3*t*gamma(5/2)*cos(x/4)*sin(x/4)^7)/(5625*pi^(1/2)) + (64*e3*t*gamma(5/2)*cos(x/4)^7*sin(x/4))/(5625*pi^(1/2)) - (64*e3*t*gamma(5/2)*cos(x/4)^3*sin(x/4)^5)/(5625*pi^(1/2)) + (64*e3*t*gamma(5/2)*cos(x/4)^5*sin(x/4)^3)/(5625*pi^(1/2))-u*uxxx-u*ux-3*ux*uxx)^2;
%ur02 = ((4*e2*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^6)/(75*gamma(2/5)) - (4*e2*t^(1/5)*gamma(7/5)/gamma(6/5)*sin(x/4)^6)/(75*gamma(2/5)) - (4*e3*t^(1/5)*gamma(7/5)/gamma(6/5)*sin(x/4)^6)/(75*gamma(2/5)) + (4*e3*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^6)/(75*gamma(2/5)) - (16*e1*gamma(6/5)*cos(x/4)*sin(x/4)^3)/(15*gamma(1/5)) - (16*e1*gamma(6/5)*cos(x/4)^3*sin(x/4))/(15*gamma(1/5)) - (16*e3*t^(2/5)*gamma(8/5)/gamma(7/5)*cos(x/4)*sin(x/4)^7)/(1125*gamma(3/5)) + (16*e3*t^(2/5)*gamma(8/5)/gamma(7/5)*cos(x/4)^7*sin(x/4))/(1125*gamma(3/5)) - (4*e2*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^2*sin(x/4)^4)/(75*gamma(2/5)) + (4*e2*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^4*sin(x/4)^2)/(75*gamma(2/5)) + (4*e3*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^2*sin(x/4)^4)/(25*gamma(2/5)) + (4*e3*t^(1/5)*gamma(7/5)/gamma(6/5)*cos(x/4)^4*sin(x/4)^2)/(15*gamma(2/5)) - (16*e3*t^(2/5)*gamma(8/5)/gamma(7/5)*cos(x/4)^3*sin(x/4)^5)/(1125*gamma(3/5)) + (16*e3*t^(2/5)*gamma(8/5)/gamma(7/5)*cos(x/4)^5*sin(x/4)^3)/(1125*gamma(3/5))-u*uxxx-u*ux-3*ux*uxx)^2;
resu = expand(ur08);
J = int(int(resu,x,0,2*pi),t,0,1)
