function qdot = quadode(t,q,T,f1,f2,d,J,m,g)
%q is the initial condition state vector : [phi theta v xposition w yposition] for ODE45

%fix(t) takes the integer value from decimal value
    omprime = (1/J)*T(fix(t));
    thprime = q(1); %from initial condition
    vprime = (1/m)*(f1(fix(t))+f2(fix(t)))*sin(q(2));
    xprime = q(3);
    wprime = (1/m)*((f1(fix(t))+f2(fix(t)))*cos(q(2))-m*g); 
    yprime = q(5);

qdot = [omprime; thprime; vprime; xprime; wprime; yprime];
%qdot = [1;1;1;1;1;1];
