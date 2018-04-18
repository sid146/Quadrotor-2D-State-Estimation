clc, clear all, close all
A = [0 0 0 1 0 0;
    0 0 0 0 1 0 ;
    0 0 0 0 0 1;
    0 0 1 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0]

B = [0 0;
    0 0;
    0 0;
    0 0;
    1 0;
    0 1;]

CTRB = [B A*B A^2*B A^3*B A^4*B A^5*B]

rank(CTRB)

rank(ctrb(A,B))

% thus Rank = 6;

%%
%finding Ko

%dx =    A x     +    B u
%       6x6 6x1      6x2  2x1

%   u=   -  k   x               (x represents output x)
%  2x1     2x6    6x1


K0 = zeros(2,6)

%%
clc, clear all, close all

n=2;
% Generate a dense n x n symmetric, positive definite matrix

A = rand(n,n); % generate a random n x n matrix

A = [.1 10;10 .1];
% construct a symmetric matrix using either
A = 0.5*(A+A'); 
% or A = A*A';
% The first is significantly faster: O(n^2) compared to O(n^3)

% since A(i,j) < 1 by construction and a symmetric diagonally dominant matrix
%   is symmetric positive definite, which can be ensured by adding nI
A = A + n*eye(n);

eig(A)

% usage of R :  u'Ru so (2x1)'R(2x1), so R=2x2

%% 3 (b) -> optimal control usin ARE

m=2; n=6; rN=75;

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
plot (Nrm(1:15))
ylabel('norm(p_i_+_1 - p_i) / norm(p_i_+_1)')
xlabel('iterations')


for i=1:length(rP)
    Z(i) = norm(rP(:,:,i));
end

%%  optimal time varying 

Kb= shiftdim(rK)
K2= shiftdim(Kb)
figure(2)
subplot(3,1,1)
plot(Kb(1,:))
ylabel('K1')
hold on;

subplot(3,1,2)
plot(Kb(2,:))
ylabel('K2')
hold on;

subplot(3,1,3)
plot(Kb(3,:))
ylabel('K3')
xlabel('iterations')
hold on;






