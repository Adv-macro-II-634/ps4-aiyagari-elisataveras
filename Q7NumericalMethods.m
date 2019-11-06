%PART 7
% USING ALTERNATIE VERSION

clear
close all

% PARAMETERS
alpha=1/3;
beta = .99; % discount factor 
sigma = 2; % coefficient of risk aversion

% ASSET
a_bar = 0; %lower bound of grid points
rho=0.5;  % productivity shock persistence
sigma_e=0.2; 


% discretizing z Rouwenhorst's method
m= 5;

[zgrid,P]=rouwenhorst(rho,sigma_e,m);
% invariant distribution 
P1=P^100; 
 Pinv = P1(1,:)';
%agregate invariate labor 
Ns=exp(zgrid*Pinv);
zgrid=exp(zgrid);

%discretizing the assets
a_lo=a_bar;
a_hi = 3;%upper bound of grid points
num_a = 500;
a = linspace(a_lo, a_hi, num_a); % asset (row) vector

% INITIAL GUESS FOR CAPITAL
K_min=0;
K_max=20;
K_guess=(K_min + K_max) / 2;

n=1;
aggsav = 1 ;

%1) solving the value function on a coarse grid, and using the result as a starting value (use linear interpolation)

%2)  policy function iteration (use k = 30 policy iterations for each opti- mization step)
% here is that I will use the Q matrix
k=30; %policy iterations for each opti- mization step
l=



%3)  linear interpolation
%4) extra credit: cubic spline interpolation.