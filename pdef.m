function pd=pdef(size)

% Generate a dense n x n symmetric, positive definite matrix

pd = rand(size,size); % generate a random n x n matrix

% construct a symmetric matrix using either
pd = 0.5*(pd+pd'); 
% or A = A*A';
% The first is significantly faster: O(n^2) compared to O(n^3)

% since A(i,j) < 1 by construction and a symmetric diagonally dominant matrix
%   is symmetric positive definite, which can be ensured by adding nI
pd = pd + size*eye(size);
pd;
