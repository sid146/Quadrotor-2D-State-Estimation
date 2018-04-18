clear all; close all; clc;
%----Exact Solution
ysol = @(t) t/(0.5 + log(t));

%----Defining f(t,y(t))
f = @(t,y) (t*y - y.^2)/t^2 ;

%----Defining Stepsize
h = 1/4;
dt=h;

%----Using Crank Nicolson Method
tcr = (1:dt:3)';
ycr = zeros(size(tcr,1),1);
ycr(1) = 2;
for ii = 2:size(tcr,1)
% Scheme is y_{n+1} = y_n +dt/2*(f(t_n,y_n)+f(t_{n+1},y_{n+1})
myequation_cr = @(x) x-ycr(ii-1) - dt/2*(f(tcr(ii-1),ycr(ii-1)) + f(tcr(ii),x));
ycr(ii) = fzero(myequation_cr,ycr(ii-1));
end

%----1) Using 3rd Order Runge-Kutta Method
t_rk3 = (1:dt:3)';
y_rk3 = zeros(size(t_rk3,1),1);
y_rk3(1) = 2;
for ii = 2:length(t_rk3)
% Scheme is y_{n+1} = y_n + 1/9*(2*k_1 + 3*k_2 + 4*k_3)
k1 = dt*f(t_rk3(ii-1),y_rk3(ii-1));
k2 = dt*f(t_rk3(ii-1)+0.5*dt,y_rk3(ii-1)+0.5*k1);
k3 = dt*f(t_rk3(ii-1)+0.75*dt,y_rk3(ii-1)+.75*k2);
y_rk3(ii) = y_rk3(ii-1) + 1/9*(2*k1 + 3*k2 + 4*k3);

end

save A11.dat y_rk3 -ascii;

%----1.2.1) Y_c
for i=1:9
    h_s(i) = 2^-(i-1);
end

%----Exact Solution
for j = 1:length(h_s)

    t_m = (1:h_s(j):3)';
    for i=1:length(t_m)
    dt = h_s(j);
    y_ex(i) = ysol(t_m(i)); %Exact solution
    y_exm = y_ex';
    end 
    
    %----Crank Nicolson
    ycr = zeros(size(t_m,1),1);
    ycr(1) = 2;
        for ii = 2:size(t_m,1)
        % Scheme is y_{n+1} = y_n +dt/2*(f(t_n,y_n)+f(t_{n+1},y_{n+1})
        myequation_cr = @(x) x-ycr(ii-1) - dt/2*(f(t_m(ii-1),ycr(ii-1)) + f(t_m(ii),x));
        ycr(ii) = fzero(myequation_cr,ycr(ii-1));
        end
    
    %----RK3
    y_rk_3 = zeros(size(t_m,1),1);
    y_rk_3(1) = 2;
        for jj = 2:size(t_m,1)
        % Scheme is y_{n+1} = y_n + 1/9*(2*k_1 + 3*k_2 + 4*k_3)
        k1 = dt*f(t_m(jj-1),y_rk_3(jj-1));
        k2 = dt*f(t_m(jj-1)+0.5*dt,y_rk_3(jj-1)+0.5*k1);
        k3 = dt*f(t_m(jj-1)+0.75*dt,y_rk_3(jj-1)+.75*k2);
        y_rk_3(jj)= y_rk_3(jj-1) + 1/9*(2*k1 + 3*k2 + 4*k3);
        
        end
    
 %----Error Matrix
 error_m(j,:) = [dt,norm(ycr-y_exm,inf),norm(y_rk_3-y_exm,inf)];

end

for i=1:8
r_trap(i) = error_m(i,2)/error_m(i+1,2);
r_rk3(i) = error_m(i,3)/error_m(i+1,3);
end

ttrap = r_trap';
trk3 = r_rk3';
save A12.dat ttrap -ascii;
save A13.dat trk3 -ascii;

%
% Exercise 2
%

m=80;
G=9.81;
tspan = 0:0.01:5;
y0 = [600,0];

%[t,y,tfinal] = ode45(@(t,y)[y(2);-g+a/m],[0,Inf],y0,opts)
[t,y]=ode45(@(t,y)paratrooper(t,y,m,G),tspan,y0);
tmp=y(length(y),1);
save A21.dat tmp -ascii;

%---- Finding Tfinal
opts = odeset('events', @g);
[t,y,tfinal]=ode45(@(t,y)paratrooper(t,y,m,G),[0, Inf],y0,opts);
save A22.dat tfinal -ascii;

%
% Exercise 3
%
%----1) Using 4th Order Runge-Kutta Method
u=5;
dt=0.125;
t_rk4 = (0:dt:50);
x_rk4 = zeros(2,length(t_rk4));
x_rk4(:,1)=[1;1];

for ii = 2:length(t_rk4)
% % Scheme is y_{n+1} = y_n + 1/6*(k_1 + 2*k_2 + 2*k_3+ k_4)
% K_1 = dt*fun(t_rk4(ii-1),x_rk4(:,ii-1));
% K_2 = dt*fun(t_rk3(ii-1)+0.5*dt,y_rk3(ii-1)+0.5*k1);
% K_3 = dt*fun(t_rk3(ii-1)+0.75*dt,y_rk3(ii-1)+.75*k2);
% myequation_rk3 = @(x) x-y_rk3(ii-1) - 1/9*(2*k1 + 3*k2 + 4*k3);

k_1 = dt*fun(t_rk4(ii-1),x_rk4(:,ii-1));
k_2 = dt*fun(t_rk4(ii-1)+0.5*dt,x_rk4(:,ii-1)+0.5*k_1);
k_3 = dt*fun(t_rk4(ii-1)+0.5*dt,x_rk4(:,ii-1)+.5*k_2);
k_4 = dt*fun(t_rk4(ii-1)+dt,x_rk4(:,ii-1)+k_3);
x_rk4(:,ii)=x_rk4(:,ii-1) + 1/6*(k_1 + 2*k_2 + 2*k_3 + k_4);

end
 
x=x_rk4(1,:);
x=x';

save A31.dat x -ascii;

