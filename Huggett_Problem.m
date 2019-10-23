% PROGRAM NAME: ps4huggett.m
clear, clc

% PARAMETERS
beta = .9932; % discount factor 
sigma = 1.5; % coefficient of risk aversion
b = 0.5; % replacement ratio (unemployment benefits)
y_s = [1, b]; % endowment in employment states
PI = [.97 .03; .5 .5]; % transition matrix


% ASSET VECTOR
a_lo = -2; %lower bound of grid points
a_hi = 5;%upper bound of grid points
num_a = 50;

a = linspace(a_lo, a_hi, num_a); % asset (row) vector

% INITIAL GUESS FOR q
q_min = 0.9;
q_max = 1;
q_guess = (q_min + q_max) / 2;

% ITERATE OVER ASSET PRICES
aggsav = 1 ;
while abs(aggsav) >= 0.01 
    q_guess = (q_min+q_max)/2;
    % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus, a', q_guess * a);
    cons = bsxfun(@plus, cons, permute(y_s, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret(cons < 0) = -inf;
    
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(2, num_a);
    
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    iteration = 0;
    while v_tol >.0001
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        v_mat = ret + beta * repmat(permute(PI * v_guess,[3 2 1]), [num_a, 1 ,1]);
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        
        [vfn, pol_indx] = max(v_mat, [], 2);
        vfn = shiftdim(vfn,2);
        % The distance between value function and v_guess
        v_tol = max(abs(vfn - v_guess));
        v_guess = vfn;
        iteration = iteration +1;
  
    end
    
    % KEEP DECSISION RULE
    %pol_fn = a(pol_indx);
    pol_indx = permute(pol_indx, [3 1 2]);
    % policy function
    pol_fn = a(pol_indx)
    
    % SET UP INITITAL DISTRIBUTION
    Mu = ones(2,num_a)/(2*num_a);
   
    
    
      %dis= 1;
      iter = 0;
     %while dis > 0.01
     for i = 1 :10000
    % ITERATE OVER DISTRIBUTIONS
    [state, a_num, mass] = find(Mu); % find non-zero indices
    
    iter = iter +1;
    MuNew = zeros(size(Mu));
    for ii = 1:length(state)
        apr_num = pol_indx(state(ii), a_num(ii)); % which a_num to go next period by policy fn?
        MuNew(:, apr_num) = MuNew(:, apr_num) + (PI(state(ii), :) * mass(ii))' ;% which mass of households goes to which exogenous state?
        dis = max(abs(Mu(:) - MuNew(:)));
        Mu = MuNew;
    end
    end
    % Find stationary distribution
    Mu;
    
     % Check Market Clearing
     % integral
     agg=Mu.* pol_fn;
     aggsav = sum(agg(:));
     if aggsav > 0.1
         q_min = q_guess;
     else
         q_max = q_guess;
     end
         
end
q_final = q_guess

% wealth distribution graph
x = 1:1:50;
plot(x,Mu(1,:),'k')
hold on
plot(x,Mu(2,:),'r')
hold off

% generate population according to Mu
pop = ones(1,10000);
% calculate for income
val_income = ones(1,10000);
val_income(1:566) = 0.5;
val_income = val_income -0.6;

[g_income,l] = mygini(pop,val_income,1);

% Gini for income
gini_income = g_income


% wealth
pop = ones(1,10000);
a_new = [a;a];
a_new(1,:) = a_new(1,:) + 1;
a_new(2,:) = a_new(2,:) + 0.5;
a_new = [a_new(1,:) , a_new(2,:)];
Mu_dist = [Mu(1,:),Mu(2,:)];
wealth_pop = Mu_dist;
wealth_val = a_new;
[g_wealth,l_wealth] = mygini(wealth_pop,wealth_val,1);

