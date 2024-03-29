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
pfi=30;
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
    
     pol_indx = permute(pol_indx, [3 1 2]);
    % KEEP DECSISION RULE        
    g = a(pol_indx); % policy function
    %what is the z value of the policy function
     [row,col] = find(g);
    %Z = zgrid(pol_indx); % policy function
    Z=repmat(zgrid,1,num_a);
    %calculating return function 
    conspf = bsxfun(@minus, r  * a', g');
    conspf = bsxfun(@plus, conspf, (Z'*w));
    retpf = (conspf .^ (1-sigma)) ./ (1 - sigma); % current period utility  
    retpfvec = retpf(:);
    retpfvec_con = v_guess(:);
    Q = makeQmatrix(pol_indx, P);
        for pii = 1:pfi
            retpf_new = retpfvec+beta*Q*retpfvec_con;
            retpfvec_con = retpf_new;
        end;
    v_guess = reshape(retpfvec_con, m, num_a);
        i=i+1;
    end;


 
   
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