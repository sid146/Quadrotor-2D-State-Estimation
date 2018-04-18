function xdot=funNL(t,x,ko)

g=9.81;m=1;I=1;w=1;
% u = U(:,t);
unl= [m*g;0]-ko*x;

h=x(1);
v=x(2) ;
th=x(3) ;
dh=x(4)  ;
dv=x(5) ;
dth=x(6);
u1=unl(1);
u2=unl(2);

dot2h = (u1*sin(th))/m;
doth  = dh;
dot2v = -m*g +(cos(th)*u1);
dotv = dv;
dot2th = (u2*1/I);
dotth = dth;
xdot = [doth;dotv;dotth;dot2h;dot2v;dot2th];
