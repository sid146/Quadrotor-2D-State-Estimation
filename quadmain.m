clc, clear all, close all

g=9.81;m=1;I=1;

 x0 = [0;0;0;0;0;0]';
% x0 = [0.1 0.1 0.1 0 0.1 0]'; 
Qxcl(:,1) = x0; 
N=10000; %number of time steps
t=linspace(0,10, N+1); %taking the time step
step = t(2)-t(1);
w=5;

Qx = zeros(6, N+1); %initial condition is set at x[:,0] to [0;0;0;0;0;0]
for i =1:N+1
   u(:,i)= [(m*g + sin(2*t(i)*pi*w));0]; 
end
% Simulation forward euler:
for i=2:N+1
    % defining the U vector
    Qx(:,i) = step*(fun(i,Qx(:,i-1),[(m*g + sin(2*t(i)*pi*w));0]))+ Qx(:,i-1);
end

% % ode 45 verification
% [tq,q]=ode45(@(tq,q)fun(tq,q),t,Qx(:,1));
% q=q';

% plot(q(2,:))
% hold on;  
% plot(q(5,:))

% figure(1)
% hold on
% plot(q(1,:))
% hold on;
% plot(q(2,:))
% hold on;
% plot(q(5,:))

figure(2)

subplot(3,1,1)
plot(Qx(1,1:5000))
ylabel('h')
hold on;
title('Forward Euler Simulation with w=5')
subplot(3,1,2)
plot(Qx(2,1:5000))
hold on;
ylabel('v')
subplot(3,1,3)
plot(Qx(3,1:5000))
ylabel('theta')
xlabel('time steps');

figure(3)
subplot(3,1,1)
plot(Qx(4,1:5000))
ylabel('dh')
hold on;
subplot(3,1,2)
plot(Qx(5,1:5000))
hold on;
ylabel('dv')
subplot(3,1,3)
plot(Qx(6,1:5000))
ylabel('dtheta')
xlabel('time steps');

figure(4)
subplot(2,1,1)
plot(u(1,1:5000))
hold on;
ylabel('u1')
subplot(2,1,2)
plot(u(2,1:5000))
ylabel('u2')
xlabel('time steps');


%% Part 2 : Stabilisation:

%Linearised Matrices

A = [0 0 0 1 0 0;
    0 0 0 0 1 0 ;
    0 0 0 0 0 1;
    0 0 9.8 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0]

B = [0 0;
    0 0;
    0 0;
    0 0;
    1 0;
    0 1;]

CTRB = [B A*B A^2*B A^3*B A^4*B A^5*B]

C = [1,0,0,0,0,0;
    0,1,0,0,0,0;
    0,0,1,0,0,0;
    0,0,0,1,0,0;
    0,0,0,0,1,0;
    0,0,0,0,0,1]



            % Ans : CTRB =
            % 
            %       [  0   0   0   0   0   0   0  9.8  0   0   0   0
            %          0   0   1   0   0   0   0   0   0   0   0   0
            %          0   0   0   1   0   0   0   0   0   0   0   0
            %          0   0   0   0   0  9.8  0   0   0   0   0   0
            %          1   0   0   0   0   0   0   0   0   0   0   0
            %          0   1   0   0   0   0   0   0   0   0   0   0  ]


rank(CTRB)
% Ans : 6; verified with rank(ctrb(A,B))

% thus Rank = 6;

% ---> 2(c)(1)
D=eig(A);
lambda = rand(10);
lambda= lambda(1,:);
for k=1:length(lambda)
rtest(k)=rank(ctrb(-lambda(k)*eye(6) - A,B));
end
% satisfies the controllability condition

% ---> 2(c)(2)

p = [-3.5; -3.6; -3.7; -3.8; -4.5; -4.9]; % be the mu values (the shift in poles)
K0 = place(A,B,p); % K0 from matlab solver


%(-muI-A)W + W(-muI-A)' + BB' =  lyapunov
        lA = (diag(p) - A);
        lQ = B*B';
        W = lyap(lA,lQ);

% from the calculation from 2(C)(3) we know that inv(P) = W; sp P = inv(W)

        lP = inv(W);

        K0= 0.5*B'*lP ;

%checking controllability
        eig(A-B*K0)
                % ans = yes controllable!
                % 
                %   -4.0094 + 8.9135i
                %   -4.0094 - 8.9135i
                %   -3.9406 + 0.8051i
                %   -3.9406 - 0.8051i
                %   -4.0500 + 3.5718i
                %   -4.0500 - 3.5718i
                
% closed Loop A
Acl = A-B*K0;
                
%% -----> 2(d) Simulation of closed loop dynamics from a small perturbation
Qxlcl = zeros(6, N+1);
Qxcl = zeros(6, N+1); %initial condition is set at x[:,0] to [0;0;0;0;0;0]
% Qxcl(:,1) = [0.1 0.1 0.1 0 0.1 0]'; 
Qxcl(:,1) = [.5 .4 .3 .3 .5 .7]';
Qxlcl(:,1) = Qxcl(:,1);
unl = zeros(2, N+1);
% Simulation forward euler:
for i=2:N+1
    % defining the U vector
    unl(:,i) = -K0*Qxcl(:,i-1)+ [9.8;0];
     Qxcl(:,i) = step*(fun(i,Qx(:,i-1),unl(:,i)))+ Qx(:,i-1);
 Qxlcl(:,i) = step*(clloop(t(i),Qxlcl(:,i-1)))+ Qxlcl(:,i-1);
end

% ode 45 verification
% [tq,qcl]=ode45(@(tq,qcl)clloop(tq,qcl),t,Qxcl(:,1));
% qcl=qcl';
% ode 45 verification
[tq,qclo]=ode45(@(tq,qclo)funNL(tq,qclo,K0),1:0.01:10,Qxcl(:,1));
qclo=qclo';



figure(45)
subplot(3,1,1)
plot(qclo(1,:))
ylabel('h')
title('implementing feedback on nonlinear quadcopter')
hold on;
subplot(3,1,2)
plot(qclo(2,:))
ylabel('v')
hold on;
subplot(3,1,3)
plot(qclo(3,:))
ylabel('theta')
hold on;

figure(46)
subplot(3,1,1)
plot(qclo(4,:))
ylabel('dh')
hold on;
subplot(3,1,2)
plot(qclo(5,:))
ylabel('dv')
hold on;
subplot(3,1,3)
plot(qclo(6,:))
hold on;
ylabel('dth')
xlabel('time step')

%input calculation
unlv=diag([m*g;0])*ones(2,length(qclo));
unlv= unlv-K0*qclo;

figure(47)
subplot(2,1,1)
plot(unlv(1,:))
ylabel('u1')
hold on;
subplot(2,1,2)
plot(unlv(2,:))
ylabel('u2')



% % 
% figure(3)
% hold on
% plot(Qxlcl(1,1:2000))
% hold on;
% plot(Qxlcl(2,1:2000))
% hold on;
% plot(Qxlcl(3,1:2000))
% hold on;
% plot(Qxlcl(4,1:2000))
% hold on;
% plot(Qxlcl(5,1:2000))
% hold on;
% plot(Qxlcl(6,1:2000))
% title('closed loop response of linearised system to small perturbation');
% xlabel('time step')
% legend ('h','v','th','dh','dv','dth')
% 
% % figure(4)
% % 
% subplot(3,1,1)
% plot(Qx(1,1:5000))
% ylabel('h')
% hold on;
% title('closed loop response of non linear system to small perturbation')
% subplot(3,1,2)
% plot(Qx(2,1:5000))
% hold on;
% ylabel('v')
% subplot(3,1,3)
% plot(Qx(3,1:5000))
% ylabel('theta')
% xlabel('time steps');
% 
% figure(5)
% subplot(3,1,1)
% plot(Qx(4,1:5000))
% ylabel('dh')
% hold on;
% subplot(3,1,2)
% plot(Qx(5,1:5000))
% hold on;
% ylabel('dv')
% subplot(3,1,3)
% plot(Qx(6,1:5000))
% ylabel('dtheta')
% xlabel('time steps');
% 
% figure(6)
% subplot(2,1,1)
% plot(u(1,1:5000))
% hold on;
% ylabel('u1')
% subplot(2,1,2)
% plot(u(2,1:5000))
% ylabel('u2')
% xlabel('time steps');
%% 3 (b) -> optimal control usin ARE

m=2; n=6; rN=N; %rN =75

rQ = pdef(n);
rR = pdef(m)


rP = zeros(n,n,rN+1); % time goes from zero to N, so indices from 1 to N+1
rP(:,:,rN+1) = rQ; % we start by setting the final P value as Q
rK = zeros(m,n,rN+1);
Nrm = zeros(rN,1);
for i = rN:-1:1

% algebraic riccatti recursion eqn for P
rP(:,:,i) = rQ + A'*rP(:,:,i+1)*A - A'*rP(:,:,i+1)*B*pinv(rR+B'*rP(:,:,i+1)*B)*B'*rP(:,:,i+1)*A;

% optimal K_t( time varying feedback) associated with each iteration
rK(:,:,i) = -pinv(rR+B'*rP(:,:,i+1)*B)*B'*rP(:,:,i+1)*A;

%norm gives a scalar value for the P matrix. (2norm is used here to show
%convergence
Nrm(rN+1-i) = norm(rP(:,:,i+1)-rP(:,:,i))/norm(rP(:,:,i+1)); 
end

%showing convergence of P
figure(9)
plot (Nrm(1:15))
ylabel('norm(p_i_+_1 - p_i) / norm(p_i_+_1)')
xlabel('iterations')


for i=1:length(rP)
    Z(i) = norm(rP(:,:,i));
end

%%  optimal time varying 
rKb = reshape(permute(rK, [2 1 3]), size(rK, 2), [])'
rk1 = zeros(76,6);
rk2 = zeros(76,6);
s=0;t=0;
for i=1:length(rKb)
    if(mod(i,2)==0)
        s=s+1;
        rk2(s,:)=rKb(i,:);
    else
        t=t+1;
        rk1(t,:)=rKb(i,:);
    end
end

        
        
figure(7)
subplot(2,1,1)
plot(rk1(:,1))
ylabel('K1')
hold on;

subplot(2,1,2)
plot(rk2(:,1))
ylabel('K2')
hold on;


%% finding initial condition x0, with norm not exceeding one, that maximises value of J


% [v,d] = eig(rP(:,:,1))

%maximizing a quadratic form subject to ?x0?=1?x0?=1 constraints for a symmetric matrix corresponds to a special matrix norm operator whose value is the maximum eigenvalue and the maximizer is the associated eigenvector.
% 
% rx0 = v(:,1);

%norm_x0 = norm(x0)
% 
% creating x vector

rx0 = [0.1 0.1 0.1 0 0.1 0]'  ;
rx = zeros(n,rN+1);
rx(:,1) = rx0;

ru = zeros(m,rN+1);
ru(:,1) = [rk1(1,:);rk2(1,:)] * rx(:,1);

% %% 2( Simulating system from x0 with u0 (normal state space simulation)
% 
for (i=1:rN)
    rx(:,i+1) = A*rx(:,i) + B*ru(:,i);
%     u(:,i+1) = Kb(:,i+1)'*x(:,i+1); % u* calculated from Kb that we calculated earlier
    ru(:,i+1) = [rk1(i+1,:);rk2(i+1,:)] * rx(:,i+1);
end

figure(3)
subplot(7,1,1)
plot(ru(1,1:60))
hold on;
title('optimal control using ARE')
plot(ru(2,1:60))
legend('u1','u2');
ylabel('u*')
hold on;

subplot(7,1,2)
plot(rx(1,1:60))
ylabel('h')
hold on;

subplot(7,1,3)
plot(rx(2,1:60))
ylabel('v')
hold on;

subplot(7,1,4)
plot(rx(3,1:60))
ylabel('theta')
hold on;

subplot(7,1,5)
plot(rx(4,1:60))
ylabel('dh')
hold on;

subplot(7,1,6)
plot(rx(5,1:60))
ylabel('dv')
hold on;

subplot(7,1,7)
plot(rx(3,1:60))
ylabel('dtheta')
xlabel('time')
hold on;

% 
%% on nonlinear quadcopter - check if required

 % rn number of time steps
% t=linspace(0,100, rN+1); %taking the time step
% step = t(2)-t(1);
w=1;
 
 
ru1 = ru+[(1*g ) 0;0 0]*ones(2,rN+1); 
 

Qx = zeros(6, rN+1); %initial condition is set at x[:,0] to [0;0;0;0;0;0]
Qx(:,1)=[0.1 0.1 0.1 0 0.1 0]';
% Simulation forward euler:
for i=2:rN+1
    
    % defining the U vector
    
     Qx(:,i) = step*(fun(i,Qx(:,i-1),ru1))+ Qx(:,i-1);
end

figure(9)
hold on
plot(Qx(1,1:80))
hold on;
plot(Qx(2,1:80))
hold on;
plot(Qx(5,1:80))


%% Part 4 Kalman Filter on closed loop system
%N=100;

%----> 4(a) Discretisation

D = zeros(6,2);
T = step; %(discretization time step)
sysc = ss(Acl,B,C,D);
sysd = c2d(sysc,T);
Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;

dA = expm(Acl*T)
dB = inv(Acl)*(dA -eye(6))*B;
dC = C;
dD = zeros(6,2)


%----> 4(b) generating noisy samples of discrete system

ksigma = .1*eye(6);
% kw = normrnd(0,0.1,N+1,1);
% kv = normrnd(0,0.2,N+1,1);

% kw1 = zeros(2,N+1);
% kv1 = zeros(6,N+1);

H = eye(6); 
F = dB;
% for i=1:N+1
%     kw1(:,i) = kw(i)*[1;1];
%     kv1(:,i) = kv(i)*[1;1;1;1;1;1];
% end

%Nonlinear Simulation
kx0 = [0.1 0.1 0.1 0 0.1 0]';


kx_real = zeros(6,N+1);
kx_real(:,1) = kx0;
for i = 2:N 
    
kx_real(:,i) = Ad*kx_real(:,i-1) ;
end


figure(11)

plot(kx_real(1,1:60))
hold on;
% plot(kx_real(2,1:60))
% hold on;
% plot(kx_real(3,1:60))
% hold on;
% plot(kx_real(4,1:60))
% hold on;
% plot(kx_real(5,1:60))
% hold on;
% plot(kx_real(6,1:60))
% ylabel('x6')
% xlabel('time')

%% simulating noisy observation

kx = zeros(6,N+1);
kx(:,1) = [0.1 0.1 0.1 0 0.1 0]';
for i = 2:N+1
    kw = normrnd(0,0.1,2,1);
    kv = normrnd(0,0.2,6,1);
kx(:,i) = Ad*kx(:,i-1) + F*kw;
kynoN(:,i) = dC*kx(:,i);
ky(:,i) = dC*kx(:,i)+ H*kv;
end
% 


plot(ky(1,1:60))
hold on;
figure(12)
subplot(6,1,1)
plot(ky(1,1:60))
ylabel('y_h')
hold on;
title('Noisy Samples of discretized system')
subplot(6,1,2)
plot(ky(2,1:60))
ylabel('y_v')
hold on;
subplot(6,1,3)
plot(ky(3,1:60))
ylabel('y_t_h_e_t_a')
hold on;
subplot(6,1,4)
plot(ky(4,1:60))
ylabel('y_d_h')
hold on;
subplot(6,1,5)
plot(ky(5,1:60))
ylabel('y_d_v')
hold on;
subplot(6,1,6)
plot(ky(6,1:60))
ylabel('y_d_t_h_e_t_a')

%% % kalman filter

Q = [1     0     0     0     0     0
     0     1     0     0     0     0
     0     0     1     0     0     0
     0     0     0     1     0     0
     0     0     0     0     1     0
     0     0     0     0     0     1];  % process noise w

R = 0.1*eye(6); % measurement noise v
K = [];


xest_b(:,1) =kx0; %before
xest_c(:,1)= kx0; %current

Pest_b(:,:,1) = ksigma;
Pest_c(:,:,1) = ksigma;

for n=1:rN

xest_b(:,n+1) = Ad*xest_c(:,n);%+dB*u; %estimate with data before current
Pest_b(:,:,n+1) = Ad*Pest_c(:,:,n)*Ad' + Q; % estimate with data before current

K(:,:,n) = Pest_c(:,:,n)*(dC'*inv(dC*Pest_c(:,:,n)*dC' + R));

xest_c(:,n+1) = xest_c(:,n)+(K(:,:,n)*(ky(:,n)-dC*xest_c(:,n)));
Pest_c(:,:,n+1) = (eye(6)- K(:,:,n)*C)*Pest_c(:,:,n);
yest(:,n) = dC*xest_c(:,n);

end
% 
% 
%  plot(1:100,ky(1,:),1:100,yest(1,:))
%    legend('true measurement','KF estimated measurement'), axis([0 100 -0.5 0.5]),xlabel('time'),ylabel('Position observations')
% 

figure(20)

plot(ky(1,1:600));
hold on
plot (yest(1,1:600))
hold on
plot(kynoN(1,1:600))
legend('y_h noisy-real','y_h -est','y_h =Cx')

figure(21)
plot(ky(2,1:600));
hold on
plot (yest(2,1:600))
hold on
plot(kynoN(2,1:600))
legend('y_v noisy-real','y_v -est','y_v =Cx')

figure(22)
plot(ky(3,1:600));
hold on
plot (yest(3,1:600))
hold on
plot(kynoN(3,1:600))
legend('y_t_h noisy-real','y_t_h -est','y_t_h =Cx')

figure(23)
plot(ky(4,1:600));
hold on
plot (yest(4,1:600))
hold on
plot(kynoN(4,1:600))
legend('y_d_h noisy-real','y_d_h -est','y_d_h =Cx')


figure(24)
plot(ky(5,1:600));
hold on
plot (yest(5,1:600))
hold on
plot(kynoN(5,1:600))
legend('y_d_v noisy-real','y_d_v -est','y_d_v =Cx')

figure(25)
plot(ky(6,1:600));
hold on
plot (yest(6,1:600))

plot(kynoN(6,1:600))
legend('y_d_t_h noisy-real','y_d_t_h -est','y_d_t_h =Cx')

%% using psd to generate new kalman estimates

Q = pdef(6);
R = pdef(6);
K = [];


xest_b(:,1) =kx0; %before
xest_c(:,1)= kx0; %current

Pest_b(:,:,1) = ksigma;
Pest_c(:,:,1) = ksigma;

for n=1:rN

xest_b(:,n+1) = Ad*xest_c(:,n);%+dB*u; %estimate with data before current
Pest_b(:,:,n+1) = Ad*Pest_c(:,:,n)*Ad' + Q; % estimate with data before current

K(:,:,n) = Pest_c(:,:,n)*(dC'*inv(dC*Pest_c(:,:,n)*dC' + R));

xest_c(:,n+1) = xest_c(:,n)+(K(:,:,n)*(ky(:,n)-dC*xest_c(:,n)));
Pest_c(:,:,n+1) = (eye(6)- K(:,:,n)*C)*Pest_c(:,:,n);
yest(:,n) = dC*xest_c(:,n);

end
% 
% 
%  plot(1:100,ky(1,:),1:100,yest(1,:))
%    legend('true measurement','KF estimated measurement'), axis([0 100 -0.5 0.5]),xlabel('time'),ylabel('Position observations')
% 

figure(26)

plot(ky(1,1:600));
hold on
plot (yest(1,1:600))
hold on
plot(kynoN(1,1:600))
legend('y_h noisy-real','y_h -est','y_h =Cx')

figure(27)
plot(ky(2,1:600));
hold on
plot (yest(2,1:600))
hold on
plot(kynoN(2,1:600))
legend('y_v noisy-real','y_v -est','y_v =Cx')

figure(28)
plot(ky(3,1:600));
hold on
plot (yest(3,1:600))
hold on
plot(kynoN(3,1:600))
legend('y_t_h noisy-real','y_t_h -est','y_t_h =Cx')

figure(29)
plot(ky(4,1:600));
hold on
plot (yest(4,1:600))
hold on
plot(kynoN(4,1:600))
legend('y_d_h noisy-real','y_d_h -est','y_d_h =Cx')


figure(30)
plot(ky(5,1:600));
hold on
plot (yest(5,1:600))
hold on
plot(kynoN(5,1:600))
legend('y_d_v noisy-real','y_d_v -est','y_d_v =Cx')

figure(31)
plot(ky(6,1:600));
hold on
plot (yest(6,1:600))
hold on
plot(kynoN(6,1:600))
legend('y_d_t_h noisy-real','y_d_t_h -est','y_d_t_h =Cx')