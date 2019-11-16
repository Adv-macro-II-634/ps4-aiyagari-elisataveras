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
delta=0.025; % depreciation rate
% discretizing z Rouwenhorst's method
m= 5;

[zgrid,P]=rouwenhorst(rho,sigma_e,m);
% invariant distribution 
P1=P^100; 
 Pinv = P1(1,:)';
%agregate invariate labor 
zgrid=exp(zgrid');
Ns=(zgrid'*Pinv);


%discretizing the assets
a_lo=a_bar;
a_hi = 90;%upper bound of grid points
num_a = 500; %500;
a = linspace(a_lo, a_hi, num_a); % asset (row) vector

% INITIAL GUESS FOR CAPITAL
K_min=20;
K_max=80;
K_guess=(K_min + K_max) / 2;


%have max iteration, because this doesn't seem to end 
maxiter=20;     
n=1;
MarkClear = 1 ;
while abs(MarkClear) >= 0.01 && n<maxiter ;
    %price associated with the guess 
    w = (1-alpha) * K_guess^(alpha)* Ns^(-alpha);
    r = (alpha)* K_guess^(alpha-1)* Ns^(1-alpha) + (1 - delta);
    
    % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus, r  * a', a);
    cons = bsxfun(@plus, cons, permute(zgrid'*w, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility  
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret(cons < 0) = -Inf;
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(m, num_a); %m state times assets states 
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    i=1;
 while v_tol >.0001;
    % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
    value_fn = ret + beta * ...
        repmat(permute((P * v_guess), [3 2 1]), [num_a 1 1]);
    
       % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
    [vfn, pol_indx] = max(value_fn, [], 2);
    vfn = permute(vfn,[3 1 2]);
    
        % what is the distance between current guess and value function
    v_tol = max(abs(v_guess(:) - vfn(:)));
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess = vfn;
    i=i+1;
 end

 pol_indx = permute(pol_indx, [3 1 2]);
    % KEEP DECSISION RULE        
    g = a(pol_indx); % policy function
    
    %get the value of asset before and after each point 
    %since my space in my matrix is as follows:0.1804
    %I use that as spacing
    g_m= a(pol_indx)-0.1804;  % index before the current one
    g_p= a(pol_indx)+0.1804;  % index above the current one
 
    a1=g(1,:);
    a2=g(2,:);
    a3=g(3,:);
    a4=g(4,:);
    a5=g(5,:);
      
      
   %use a for loop to go over the complete 
   for i = 1:length(a1)
         
    %need to do linear approvimation over the possible values 
    %first, the initial points are two
    x_low=[g_m g];
    x_high=[g g_p];
    
    %then the value function is:
    vf=vfn';
    %create 10 numbers between the two possible numbes 
    xi_low = g_m:1/10:g; 
    xi_high = g:1/10:g_p; 
       %interpolate
       interp_l = interp1(x_low,vf,xi_high); 
 end
   
        %Golden section search
    
    % ------------------------GOLDEN SECTION METHOD----------------------------
% -------------------------------------------------------------------------
% Copyright (c) 2009, Katarzyna Zarnowiec, all rights reserved 
% mailto: katarzyna.zarnowiec@gmail.com
% -------------------------------------------------------------------------

a=g_m1;                            % start of interval
b=g_p1;                            % end of interval
epsilon=0.000001;               % accuracy value
iter= 50;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k=0;                            % number of iterations


x1=a+(1-tau)*(b-a);             % computing x values
x2=a+tau*(b-a);

f_x1=f(x1);                     % computing values in x points
f_x2=f(x2);

plot(x1,f_x1,'rx')              % plotting x
plot(x2,f_x2,'rx')

while ((abs(b-a)>epsilon) && (k<iter))
    k=k+1;
    if(f_x1<f_x2)
        b=x2;
        x2=x1;
        x1=a+(1-tau)*(b-a);
        
        f_x1=f(x1);
        f_x2=f(x2);
        
        plot(x1,f_x1,'rx');
    else
        a=x1;
        x1=x2;
        x2=a+tau*(b-a);
        
        f_x1=f(x1);
        f_x2=f(x2);
        
        plot(x2,f_x2,'rx')
    end
    
    k=k+1;
end


% chooses minimum point
if(f_x1<f_x2)
    Agold=x1;
else
     Agold=x2;
end
 
   
% SET UP INITITAL DISTRIBUTION
Mu = zeros(m,num_a);
Mu(1, 4) = 1; % initial guess: everyone employed, 0 assets

% ITERATE OVER DISTRIBUTIONS
% way 1: loop over all non-zeros states
mu_tol = 1;
while mu_tol > 1e-08
    [emp_ind, a_ind] = find(Mu > 0); % find non-zero indices
    
    MuNew = zeros(size(Mu));
    for ii = 1:length(emp_ind)
        apr_ind = pol_indx(emp_ind(ii), a_ind(ii)); 
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ...
            (P(emp_ind(ii), :) * Mu(emp_ind(ii), a_ind(ii)) )';
    end
    mu_tol = max(abs(MuNew(:) - Mu(:)));
    Mu = MuNew ;
end

    aggK= sum( g(:) .* Mu(:) ); 
    
    % now I need to update the information   and check for market clearing  

    MarkClear=aggK-K_guess;

if MarkClear > 0 ;
    K_min = K_guess ;
end ;
if MarkClear < 0;
    K_max = K_guess ;
end ;

display (['r = ', num2str(r)])
display (['Aggregate desired capital = ', num2str(K_guess)]);
display (['New K_min is ', num2str(K_min), ', new K_max is ', num2str(K_max)]);
display (['New K_guess is ', num2str((K_min + K_max)/2)]);

K_guess = (K_min + K_max)/2 ;
display (' ') ;
    n=n+1;
end